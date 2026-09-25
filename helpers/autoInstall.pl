#!/usr/bin/env perl
# autoInstaller for lotus
# Copyright (C) 2014  Falk Hildebrand, Joachim Fritscher

### use "perl autoInstall.pl -condaDBinstall -lambdaIndex" to install LotuS3 on Galaxy server

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

# contact
# ------
# Falk.Hildebrand [at] gmail.com
# 

use strict;
use warnings;
use Getopt::Long qw( GetOptions );
use Cwd 'abs_path';
use File::Copy qw(move copy);
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempdir);
use POSIX qw(uname);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3 qw(open3);
use File::Glob qw(bsd_glob); #unlike glob(), does not split patterns at whitespace in install paths
sub addInfoLtS;sub finishAI;
#subroutines to download various DBs..
sub getGG2; sub getSLV;sub getHITdb; sub getPR2db;sub getKSGP;sub getbeetax;
sub buildIndex;
sub get_DBs;
sub getS2;
sub getInfoLtS;
sub getInstallVer;
sub compile_sdm;
sub compile_LCA;
sub install_bundled_rtk;
sub version_is_newer;
sub check_version;
sub user_options;
sub command_exists;
sub run_cmd;
sub capture_cmd;
sub read_user_input;
sub gunzip_file;
sub write_config_atomic;
sub copy_file_atomic;
sub replace_tree_atomic;
sub ensure_dir;
sub verify_sha256;
# Only core Perl modules are used. Downloads go through wget or curl, or HTTP::Tiny (core)
# when the optional SSL modules it needs for https are installed.
BEGIN { die "The LotuS3 installer needs Perl 5.14 or newer (found $]).\n" if $] < 5.014 }
require HTTP::Tiny;
#can_ssl only exists in newer HTTP::Tiny releases (Perl >= 5.22)
sub http_tiny_https { return (HTTP::Tiny->can('can_ssl') && HTTP::Tiny->can_ssl) ? 1 : 0 }
my $forceUpdate=0;
my $condaDBinstall=0;
my $downloadLmbdIdx = 0; #download lambda index from webpage
my $compile_lambda=0;
my $usearchInstall = "";
my $noTelemetry = 0;
my $ontOnly = 0;
my $withBarbell = 0; #accepted for compatibility; Barbell is included with ONT tools
my $showHelp = 0;
my @pendingLambdaIndex; #databases waiting for lambda3, which a fresh install adds after the databases

# SHA-256 of every file the installer downloads, recorded from the upstream files on 2026-09-25.
# getS2 refuses any URL that is neither listed here nor given an explicit checksum by its caller,
# so a changed or substituted file stops the installation instead of being installed.
# To move to a new release: download it, check it, and replace URL and checksum together.
my %PINNED_SHA256 = (
	# programs
	'https://lotus2.earlham.ac.uk/lotus/packs/ITSx_1.1.4.tar.gz' => '513a83a10f991f83571f466e9b81b7f4fcb5958396d0cdfe9977e24ebc330ac2',
	'https://lotus2.earlham.ac.uk/lotus/packs/ncbi-blast-2.2.29+-x64-linux.tar.gz' => 'f53e97aaff424c2583cebff76a318f36610b5444657cee2b2072513ca907fedf',
	'https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Linux-x86_64.tar.xz' => '854e42d41521c483a0e6a92730fc1daf058ca8df9c0fe7af40b66f4e885d2d49',
	'https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Darwin-x86_64.zip' => '51f546f78300277962af810e66ee8019a232d99c6dd5b2e265d02c39c24241b0',
	'https://lotus2.earlham.ac.uk/lotus/packs/swarm2.1.13.zip' => 'bbfb326adc68d1ac96603a8eb828dc591dac9f1410cbac144054f078d8a433d9',
	'https://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-macosx-intel.tar.gz' => '87adbf63c61f66127f14823f1a17451fa960567d0dfdc1d72735a16d41a6a171',
	'https://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz' => 'ad062059dfff6450d4a9417846e2695087f5b7c13d64a6deef4a5bcb51b75f3e',
	'https://lotus2.earlham.ac.uk/lotus/packs/VXtractor/vxtractor.pl' => '76f7d103470a89d4dfad19a4c444b7095feaba37b7de241bc2619fc849f1f739',
	'https://lotus2.earlham.ac.uk/lotus/packs/VXtractor/HMMs.zip' => 'eb5e2dac2b0919251897510dd88405a0ac89cff4ee3457c9d581ccb08556bfe8',
	'https://github.com/iqtree/iqtree2/releases/download/v2.1.1/iqtree-2.1.1-MacOSX.zip' => '300a696f2527cadc05f87c4cfea9b2a29964b6ae754bfc52e49375f9fa241787',
	'https://github.com/iqtree/iqtree2/releases/download/v2.1.1/iqtree-2.1.1-Linux.tar.gz' => '594f23ee2ec04bfb7126c0d95f9c75efcda4718816d2a7f113689bf617fa7313',
	'https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip' => 'c388f9bb85b8ccbdb7d425ac01b5865c84ce27f2051e86fad5f9eb118dd31b1a',
	'https://mafft.cbrc.jp/alignment/software/mafft-7.471-linux.tgz' => 'f85019489117dc554da5cf43ec331ae61dc93be57d43ec35bbeaebc4f55b4380',
	'https://lotus2.earlham.ac.uk/lotus/packs/FastTree.c' => 'da148297bb64711e43e38481186228496d418bb4ec0166e09df62a72248085a0',
	'https://lotus2.earlham.ac.uk/lotus/packs/cd-hit_git.zip' => 'ef5ccb22b0d0faa816f744dbb83c5571db0ac089e49f0a6f57848ada37ac29b3',
	'https://lotus2.earlham.ac.uk/lotus/packs/rdp_classifier_2.12.zip' => '977896248189a1ce2146dd3a61d203c3c6bc9aa3982c60332d463832922f7d0a',
	'https://lotus2.earlham.ac.uk/lotus/packs/clustalo-1.2.0-Ubuntu-x86_64' => '2b04eef987d1c5ae73fafc1bd3998250607030eb7e936075d7e4edc9992e228c',
	'https://github.com/rcedgar/usearch12/releases/download/v12.0-beta1/usearch_linux_x86_12.0-beta' => '4193abead8c7e1609dd28148bb36ad9667c67647c6f784f2bdd72af9de27f3dc',
	'https://github.com/torognes/vsearch/releases/download/v2.32.0/vsearch-2.32.0-linux-x86_64.tar.gz' => 'c9d7ad4e10e942286ad84004a913e0f2f82957f6156c97a7d6d2389750dd6e41',
	'https://github.com/torognes/vsearch/releases/download/v2.32.0/vsearch-2.32.0-linux-aarch64-static.tar.gz' => 'b651509e69b5c0667cebb8e9f555bf18443e48e9decf84a3467e41fb091d74cb',
	'https://github.com/torognes/vsearch/releases/download/v2.32.0/vsearch-2.32.0-macos-universal.tar.gz' => 'a55e5c2a9a66e8dbc58543a2afe161631e0424b6a6c323eda2f87269c565aed3',
	# reference and chimera databases
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/phiX.fasta' => '398563e14ebb13248eaaa3bdcb95a36a90172e016f7c1c0527382a590d2c811a',
	'https://lotus2.earlham.ac.uk/lotus/packs/rdp_gold.fa.gz' => '1b64cfb56efa27fb28325522dd2d29f0a900888b8f1a64741c5e943b71debf69',
	'https://lotus2.earlham.ac.uk/lotus/packs/SILVA_119_LSU_93.ref.fasta.gz' => '53f3d0728d8a7e4760b3170e84e790815d1f0f5874ff6e705b96a777f2e13de6',
	'https://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_sequences.fna' => '4c4cc7c1316928308c5cca179faad02e3f97ac5c6e31393e27d721161cd1e7d5',
	'https://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_taxonomy_qiime.txt' => 'dda14eb2315fbf58130cfc75d0f416cf9d362727bad40e6204c0b66992ce47ae',
	'https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/HITdb_sequences.fasta.lba.gz' => '5fc0b88214772d0e61890f5017deaab4396c0a23050722aa13136a6fd48130cb',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.fna' => '00e6dde860b5a840b904e48f02165208782aaf20679895b4183b5cff5993e80a',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.txt' => '01b0fa5c7a23450d1f3e8ed351d23b58541cfd0e64374d33fdd7964b1f2b1d49',
	'https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/beeTax.fasta.lba.gz' => '37e2a1db53559c97e7436f775953c322224a1c56849558c0e707357c6ec5fb4b',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/PR2//pr2_version_5.0.0_SSU_UTAX.fasta.gz' => '4239d2d441f8ac8e2bb6c357a425d62a8b431df1c8f4e66105ee6d76e73f9c48',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/uchime/uchime_UNITE_16_10_22.zip' => 'd5f995e3e87074d78296538140ec2f9d6ece9bac3de26ca48b818f1215c506e3',
	'https://drive5.com/utax/data/utax_rdp_16s_tainset15.tar.gz' => 'f6655af6f78e1c612d4c9c0cd750821a7ac0f590022cf7087ef1849e5df1a05f',
	'https://drive5.com/utax/data/utax_unite_v7.tar.gz' => '6c4dca38eae84b0441a8512151b9047b4d3b355270b7e4cc2b7cd7a8d84ba8d4',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_v4.0.fasta.gz' => '8dab2647fa53f83427b95e947b115442f4d9833484ebf4dd355c33f78ec1f62f',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_plus2.tax.gz' => '6a064388f0beac30a115f9c81e9b4ad2c826d3a12e541562aca4376c24a2c8cd',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGPv4.0.fasta.lba.gz' => '33e710dc9f0c72cc761eedb7fcad6e031d5e1b84fc18fe2fc85e8befb98ffa51',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_general_release_dynamic_s_all_19.02.2025.fasta.gz' => 'a51cf593618534ee642ebcc320d8ea412c645983a2c026d887a82e27c96460db',
	'https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/Lambda3/sh_refs_v10_19.02.2025.fasta.lba.gz' => '9c279ab807e79d43dc7b788cd44736671eefc5f86a5e4fc84cd86e0c698b0403',
	'https://ksgp.earlham.ac.uk/downloads/greengenes2/GG2.2022.10.fasta.gz' => '38c581c7d18360aadb7504aec10b012fcbd1fb3b28ddade02c925a3cf0f0d6bf',
	'https://ksgp.earlham.ac.uk/downloads/greengenes2/GG2.2022.10.tax.gz' => '673122791d3f6fce079fae1c775f670775437bfa9658238278305c6399b6fc91',
	'https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz' => '7078a4e54ee962ca3108e776b8935795c744869e73e24c4718b8b5ff240410d0',
	'https://ftp.arb-silva.de/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz' => '887627935d83cc0ad88fb3cac93b154470304fec6e73cbb2a04ae625427a203d',
	'https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_LSURef_tax_silva.fasta.gz' => '43e3b8183c11343df2b0fcb98cdae5aa9bdc98f1091974698698969e5bb6a748',
	'https://ftp.arb-silva.de/release_138.1/Exports/taxonomy/tax_slv_lsu_138.1.txt.gz' => 'bfc32b01b1285857bdbc0d06bd0f62e4d8efc80c71adcad0a9c07aab0bc158e6',
	'https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/GG2.2022.10.fasta.lba.gz' => '764dfa7a04067ef6a30a92d7689dc1bec1a2256ff3b46f6dc9c72ec705f700f5',
	'https://ksgp.earlham.ac.uk/lambdaDBs/v3.0/SLV_138.1_SSU.fasta.lba.gz' => 'c320cfd668a6da868f0f9283f1f051c1a90d8041fa77e4fd5c60fed0609b8fce',
);

GetOptions(
	"h|help"          => \$showHelp,
	"forceUpdate"     => \$forceUpdate,
	"condaDBinstall"  => \$condaDBinstall,
	"downloadLmbdIdx" => \$downloadLmbdIdx,
	"lambdaIndex"     => \$compile_lambda,
	"link_usearch=s"  => \$usearchInstall,
	"no-telemetry"    => \$noTelemetry,
	"ont-only"        => \$ontOnly,
	"with-barbell"    => \$withBarbell,
) or die "Invalid command line options (see perl helpers/autoInstall.pl --help)\n";

if ($showHelp) {
	print <<'HELP';
LotuS3 autoinstaller: installs programs and reference databases into this LotuS3 directory
and registers them in lOTUs.cfg. Every download is checked against a pinned SHA-256 checksum.

Usage: perl helpers/autoInstall.pl [options]

  (no option)        interactive installation; a rerun offers to refresh databases/programs
  --ont-only         install or register only the ONT tools (minimap2, Savont, Barbell)
  -condaDBinstall    non-interactive download of the standard database set (Bioconda installs)
  -downloadLmbdIdx   download prebuilt Lambda indices instead of building them
  -lambdaIndex       build Lambda indices for the installed databases
  -link_usearch PATH register an existing USEARCH binary and exit
  --no-telemetry     do not send the anonymous installation ping (install ID and versions)
  --with-barbell     accepted for compatibility; Barbell is part of the ONT tools
  -h, --help         show this help

Update LotuS3 itself with "git pull" (GitHub checkout) or "conda update lotus3" (Bioconda).
HELP
	exit(0);
}

if ($forceUpdate) {
	#the old online updater fetched and unpacked an unverifiable archive over plain HTTP
	die "-forceUpdate is no longer supported: the online updater could not verify what it downloaded.\n"
		. "Update LotuS3 with \"git pull\" in a GitHub checkout or \"conda update lotus3\" for Bioconda,\n"
		. "then rerun perl helpers/autoInstall.pl to refresh programs and databases.\n";
}

if ($ontOnly && ($condaDBinstall || $compile_lambda || $downloadLmbdIdx || $usearchInstall ne "")) {
	die "--ont-only cannot be combined with database, update, or USEARCH-link modes.\n";
}

if ($withBarbell && ($condaDBinstall || $compile_lambda || $downloadLmbdIdx || $usearchInstall ne "")) {
    die "--with-barbell requires a program installation (normal install or --ont-only).\n";
}

if ($compile_lambda && $downloadLmbdIdx){
	die "Can't use both -lambdaIndex and -downloadLmbdIdx arguments together\nAborting..\n";
}

my $WGETpres = command_exists("wget") ? 1 : 0;
my $CURLpres = command_exists("curl") ? 1 : 0;

my $isMac = 0;
if ($^O eq "darwin"){
	$isMac = 1;
	print "Detected MAC.. will install LotuS3 for MAC\n";
} elsif ($^O !~ m/linux/){
	die "Unsupported operating system '$^O'. The LotuS3 installer supports Linux and macOS only.\n";
}

my $ldir = abs_path($0);
$ldir =~ s/\/[^\/]*$/\//;
if (! -e "$ldir/lotus3" && -e "$ldir/../lotus3"){#autoInstall might be in helpers/
	$ldir =~ s/\/[^\/]+\/$/\//;
}

print "installing into $ldir\n";
#die "\n\n$ldir\n\n";

#die ($ldir."\n");
my $bdir = $ldir."/bin/";
my $ddir = $ldir."/DB/";
my $finalWarning="";
my $configBackupWritten = 0;
my $onlyDbinstall = 0;
#options on programs to install..
my $installBlast = 2; my @refDBinstall = (0) x 10; my $ITSready = 1;my $getUTAX=1;$refDBinstall[8]=1;
my $installONT = 1; #ONT-only mode includes the whole set; detailed setup may opt out

#DEBUG
#get_programs();die;
#@refDBinstall = 0 x 10; $ITSready = 1;$getUTAX = 0; get_DBs();die;


#autoinstaller, test if install was done before
my @txt; my $mainCfg = "$ldir/lOTUs.cfg";
my $defCfg = "$ldir/configs/LotuS.cfg.def";
my $configReadPath = $mainCfg;
if (!-e $mainCfg ){
	die "Something wrong: can't find configs: $mainCfg and $defCfg" unless (-e $defCfg);
	$configReadPath = $defCfg;
}
open I,"<",$configReadPath or die "Cannot read $configReadPath: $!\n";
while (my $line = <I>){	push(@txt,$line);}
close I;
my $exe = ""; my $callret;
#print "$ldir/lOTUs.cfg";
#hand-edited configurations may lack these entries (e.g. for --ont-only)
my $UID = getInfoLtS("UID",\@txt,"??");
my $uspath = getInfoLtS("usearch",\@txt,"");

#usearch binary linking is handled by GetOptions above


##### TESTING / DEBUG ##########
#		@txt = getKSGP(\@txt);die;#@txt = getGG2(\@txt); 

#DEBUG
#@txt = getPR2db(\@txt);;exit;

# ONT-only installation leaves databases, R packages, and other tools untouched.
if ($ontOnly) {
    ensure_dir($bdir);
    check_ont_build_requirements();
    install_ont_programs();
    finishAI("none");
    print "Installed ONT programs and registered their paths in $mainCfg\n";
    exit(0);
}

###### GET USER OPTIONS ###################

my ($lver,$sver) = getInstallVer("");
if ($condaDBinstall){
	print "\n\nConda LotuS install: downloading the standard database set for LotuS3\n\n";
} else {
	print "\n\t####################################\n\t LotuS $lver Auto Installer script.\n\t####################################\n\n";
}
user_options();

###### END GET USER OPTIONS ###################


#prepare dirs
#system("rm -rf $bdir");
ensure_dir($bdir);
#system("rm -rf $ddir");
ensure_dir($ddir);
($lver,$sver) = getInstallVer("$ldir/sdm_src");



if ($UID eq "??"){
	$UID=int(rand(999999999));
	@txt = addInfoLtS("UID",$UID,\@txt,0);
}

# Validate full-install prerequisites before starting potentially large database
# downloads. Database-only modes deliberately skip this check.
my $install_dada = 1;
my $rscriptBin;
if (!$condaDBinstall && !$onlyDbinstall){
	$rscriptBin = command_exists('Rscript');
	if ($rscriptBin) {
		my $v = check_version($rscriptBin);
		if ($v < 4) {
			print("$0 requires Rscript version 4.0.0 or newer for dada2 and phyloseq.\nThe found version is older than 4.\nType\n  'c' to install LotuS3 without dada2 and phyloseq\n  't' to continue and try to install phyloseq only\n  'a' to abort installation process\n");
			my $instr = lc(read_user_input("the R version choice"));
			if ($instr eq "c") {
				$install_dada = 0;
			} elsif ($instr eq "t") {
				$install_dada = 1;
			} else {
				die "Installation aborted by user.\n";
			}
		}
	} else {
		die "$0 requires Rscript (version 4.0.0 or newer) for a full install. No Rscript was detected.\n";
	}
}


check_ont_build_requirements() unless $condaDBinstall || $onlyDbinstall;
#fail before gigabytes of databases are downloaded, not after
check_full_install_requirements() unless $condaDBinstall || $onlyDbinstall;

###################   database downloads ... #########################
get_DBs();

if ($condaDBinstall){
	build_pending_lambda_indexes(1);
	finishAI("d");
	print "Finished LotuS3 database install (Conda autoinstall)\nEnjoy LotuS3!\n";
	exit(0);
}

if ($onlyDbinstall){
	build_pending_lambda_indexes(1);
	finishAI("d");
	print "\n\nInstalled databases\nExiting autoinstaller..\n";
	exit(0);
}

######## get programs ####################

print "Several software packages have to be downloaded and this can take some time. Please be patient & grab a tea.\n\n";

# USEARCH placeholder for configurations that do not already have a valid path
if (!-e $uspath){
	$exe = $bdir."usearch_bin";
	@txt = addInfoLtS("usearch",$exe,\@txt,0);
}

###################   R packages ... #########################

if ($install_dada) {
	print("Install dada2 and other R packages\n");
	my $rscript = $ldir . "/helpers/autoInstall.R";
	die "Cannot find R package installer $rscript\n" unless (-f $rscript);
	my ($r_output,$r_status) = capture_cmd($rscriptBin, $rscript);
	print($r_output);
	if ($r_output =~ m/(Package .* could not be installed\. Please install it manually in your R environment\.)/){
		$finalWarning .= "$1\n";
	}
	if ($r_status != 0){
		$finalWarning .= "R package installation failed with exit status $r_status. Review the R output and install the missing packages manually.\n";
	}
}


#only binary installs after this point
#compile sdm/LCA and validate the bundled rtk binary
my $nsdmp = compile_sdm("$ldir/sdm_src");
@txt = addInfoLtS("sdm",$nsdmp,\@txt,1);
$nsdmp = compile_LCA("$ldir/LCA_src");
@txt = addInfoLtS("LCA",$nsdmp,\@txt,1);

#rtk (rarefaction) is a standalone helper that lotus3 itself does not call: a failure only warns
my $rtkPath = eval { install_bundled_rtk("$ldir/bin/rtk") };
if (defined $rtkPath) {
	@txt = addInfoLtS("rtk",$rtkPath,\@txt,1);
} else {
	my $msg = "rtk was not installed: $@";
	print $msg; $finalWarning .= $msg;
}

#download and install the remaining programs exactly once
get_programs();
build_pending_lambda_indexes(1); #lambda3 is installed by now

finishAI("");


print "\n\nInstallation script finished.\nPlease read the README for examples and references to proprietary software used in this pipeline.\n";


#After install on your system, open\n   ".$ldir."lOTUs.cfg\nand search for the entry \"usearch {XX}\".\nReplace {XX} with the absolute path to your usearch install, e.g. /User/Thomas/bin/usearch/usearch7.0.1001_i86linux32\n LotuS is ready to run.\n";

sub finishAI($){
	my ($vTag) = @_;
	#write new cfg file
	write_config_atomic("$ldir/lOTUs.cfg", \@txt);
	return if ($vTag eq "none");
	if (!$noTelemetry){
		#one short, best-effort request; its answer is ignored and a failure never stops the installer
		my $ping = "https://lotus2.earlham.ac.uk/lotus/in.php?ID=$UID&VERSION=$vTag$lver&SDMV=$sver";
		if ($WGETpres) { system("wget", "-q", "-T", "10", "-t", "1", "-O", File::Spec->devnull(), $ping); }
		elsif ($CURLpres) { system("curl", "-s", "-m", "10", "-o", File::Spec->devnull(), $ping); }
		elsif (http_tiny_https()) { HTTP::Tiny->new(timeout => 10)->get($ping); }
	}
	if ($finalWarning ne ""){
		print "################################\nWarnings occured during LotuS installation:\n".$finalWarning."\n################################\n";
	}
}
sub getInstallVer($){
	my ($sdmsrc) = @_;
	my $lver=0.1;
	run_cmd("chmod", "+x", "$ldir/lotus3") if -e "$ldir/lotus3";
	open Q,"<","$ldir/./lotus3" or die("Can't find LotuS main script file (lotus3)\n");
	while(<Q>){if (m/my.*selfID\s*=\s*\"LotuS\s(.*)\".*/){$lver=$1;last;}}
	close Q;
	my $sver=1.5;
	if ($sdmsrc ne ""){
		#the installed sdm reports its own version ("sdm 3.53 beta"); the source tree is only a fallback
		my $sdmBin = "$bdir/sdm";
		if (-x $sdmBin){
			my ($sdmV, $status) = capture_cmd($sdmBin, "-v");
			$sver = $1 if ($status == 0 && $sdmV =~ m/sdm\s+(\d+(?:\.\d+)+)/);
		}
		my $sdmF = "$sdmsrc/IO.h";
		if ($sver == 1.5 && -e $sdmF){
			open Q,"<",$sdmF or die("Can't open sdm file $sdmF\n");
			#static const float sdm_version = 0.71f;
			while(<Q>){if (m/static\s+const\s+float\s+sdm_version\s*=\s*(.*)f;/){$sver=$1;last;}}
			close Q;
		}
	}
	return ($lver,$sver);
}

sub addInfoLtS($ $ $ $){
	my ($cmd,$ex,$aref,$reqF) = @_;
	print "Installing $cmd:\n$ex\n";
	if ($reqF ==1 && ! -f $ex){print "Can't find required file $ex\nPlease check if the package was correctly downloaded.\nAborting..\n"; exit(5);}
	if ($reqF ==2 && ! -d $ex){print "Can't find required directory $ex\nPlease check if the package was correctly downloaded.\nAborting..\n"; exit(5);}
	my @txt = @{$aref};
	my $ss = quotemeta $cmd;
	#print "$ss\nXX\n";
	my $i=0; my $tagset=0;
	while ($i < @txt && $txt[$i] !~ m/^$ss\s/){
		#print $txt[$i]."\n";
		$i++;
		if ($i >= @txt){
			#die ("Could not find the entry \"$cmd\" in lotus configuration file. Aborting Installer..\n")
			print "Could not find the entry \"$cmd\" in lotus configuration file. Inserting anew..\n";
			push(@txt,""); last;
		}
	}
	$txt[$i-1] .= "\n" if $i > 0 && $txt[$i-1] !~ /\n$/;
	$txt[$i] = $cmd." ".$ex."\n";
	$i++;
	while ($i<@txt){ if ($txt[$i] =~ m/^$ss\s/){splice(@txt,$i,1) ; $i--;} $i++; last if ($i>=@txt); }
	print "done.\n";
	
	write_config_atomic("$ldir/lOTUs.cfg", \@txt);

	#DEBUG
	#print $txt[$i]."\n";
	return @txt;
}
               
sub getInfoLtS($ $;$){
	my ($cmd,$aref,$default) = @_;
	my $ss = quotemeta $cmd;
	foreach my $line (@{$aref}){
		chomp(my $copy = $line);
		return $1 if ($copy =~ m/^$ss\s+(.*)$/);
		return "??" if ($copy =~ m/^$ss\s*$/);
	}
	return $default if @_ > 2;
	die ("Could not find the entry \"$cmd\" in lotus configuration file. Aborting Installer..\n");
}
sub parse_hitdb($ $){
	my ($Dpre,$Dn) = @_;
	my @tdesign = (" k__"," p__"," c__"," o__"," f__"," g__"," s__");
	my $tmp = "$Dn.tmp.$$";
	open I,"<",$Dpre or die "Cannot open HITdb taxonomy input $Dpre: $!\n";
	open O,">",$tmp or die "Cannot write HITdb taxonomy output $tmp: $!\n";
	while (my $l = <I>){
		chomp $l;
		my @spl = split /\t/,$l;
		#print $spl[1]."\n";
		my @spl2 = split /;/,$spl[1];
		my $nline = "";
		if ($spl2[0] =~ m/Euryarchaeota|Crenarchaeota/){
			$nline = $spl[0]."\tk__Archaea;";
		} else {
			$nline = $spl[0]."\tk__Bacteria;";
		}
		for (my $i=1;$i<@tdesign;$i++){
			
			if (@spl2 >= $i && $spl2[$i-1] ne ""){ 
				my $tag = $spl2[$i-1]; chomp $tag;
				$nline .= $tdesign[$i].$tag;
			} else {
				$nline .= $tdesign[$i]."?";
			}
			$nline .=";" unless ($i == (@tdesign-1));
		}
		print O $nline."\n";
	}
	close I or die "Cannot close $Dpre: $!\n";
	close O or die "Cannot close $tmp: $!\n";
	rename($tmp, $Dn) or die "Cannot replace $Dn with $tmp: $!\n";
}

sub parse_PR2($ $){
	my ($DBin, $tout) = @_;
	
	print "Rewriting PR2 database..\n";
	
	my $taxTmp = "$tout.tmp.$$";
	open T,">",$taxTmp or die "Can't open PR2 taxonomy output $taxTmp: $!\n";
	open I,"<",$DBin or die "Can't open PR2 fasta $DBin: $!\n";
	open F,">","$DBin.tmp" or die "Can't open PR2 fasta tmp $DBin.tmp: $!\n";
	#>AB353770.1.1740_U;tax=k:Eukaryota,d:TSAR,p:Alveolata-Dinoflagellata,c:Dinophyceae,o:Peridiniales,f:Kryptoperidiniaceae,g:Unruhdinium,s:Unruhdinium_kevei
	my $noTax = 0; my $firstNoTax = "";

	while (my $l = <I>){
		chomp $l;
		if ($l =~ m/^>/){
			 my @spl = split /;tax=/,$l;
			$spl[0] =~ s/^>//;
			print F ">".$spl[0]."\n";
			if (!defined($spl[1]) || $spl[1] eq ""){ #no ";tax=": all ranks unknown
				$firstNoTax = $spl[0] if $noTax++ == 0;
				$spl[1] = "";
			}
			my $taxS = $spl[1];my $taxO="";
			foreach my $lvl ( ("k","p","c","o","f","g","s") ){
				my $taxL = "?";
				if ($taxS =~ m/$lvl:([^,]+)/){$taxL = $1;}
					
				if ($lvl ne "s"){
					$taxO .= "${lvl}__$taxL;";
				} else {
					$taxO .= "${lvl}__$taxL";
				}
			
			}
			#print "$taxS\n$taxO\n";
			#print T $spl[0]."\td__".$spl[1].";p__".$spl[2].";c__".$spl[4].";o__".$spl[5].";f__".$spl[6].";g__".$spl[7].";s__".$spl[8]."\n";
			print T "$spl[0]\t$taxO\n";
		} else {
			$l =~ s/U/T/g;
			$l =~ s/u/t/g;
			$l =~ s/[^ACTGactg]/N/g;
			print F $l."\n";
		}
	}
	close T or die "Cannot close $taxTmp: $!\n"; close I; close F or die "Cannot close $DBin.tmp: $!\n";
	if ($noTax){
		my $msg = "PR2: $noTax sequence headers have no ';tax=' annotation (first: $firstNoTax); their taxonomy is recorded as unknown.\n";
		print $msg; $finalWarning .= $msg;
	}
	rename($taxTmp, $tout) or die "Cannot replace $tout with $taxTmp: $!\n";
	unlink($DBin) if -e $DBin; move("$DBin.tmp", $DBin) or die "Cannot replace $DBin with $DBin.tmp: $!\n";
}


sub configured_lambda3 {
	my $lambdaIdxBin = "";#find where lambda is installed in
	foreach my $line (@txt){
		if ($line =~ m/^lambda3\s+(\S+)/ ) {$lambdaIdxBin = $1;}
	}
	if ($lambdaIdxBin ne "" && !-x $lambdaIdxBin){
		$lambdaIdxBin = command_exists($lambdaIdxBin) // "";
	}
	return ($lambdaIdxBin ne "" && -x $lambdaIdxBin) ? $lambdaIdxBin : "";
}

sub buildIndex($){
	my ($DBfna) = @_;
	return unless ($compile_lambda);
	die "Cannot build Lambda index: database file $DBfna is missing or empty\n" unless (-s $DBfna);
	if (configured_lambda3() eq "") {
		print "lambda3 is not installed yet; the Lambda index for $DBfna is built once it is.\n";
		push @pendingLambdaIndex, $DBfna;
		return;
	}
	build_lambda_index($DBfna);
}

#build the deferred indices; $final: lambda3 will not be installed later in this run
sub build_pending_lambda_indexes {
	my ($final) = @_;
	return unless @pendingLambdaIndex;
	if (configured_lambda3() eq "") {
		return unless $final;
		die "Cannot build Lambda indices (-lambdaIndex): lambda3 is not configured or executable.\n"
			. "Install Lambda 3 (or put lambda3 on PATH) and rerun with -lambdaIndex. Databases waiting: @pendingLambdaIndex\n";
	}
	build_lambda_index(shift @pendingLambdaIndex) while @pendingLambdaIndex;
}

sub build_lambda_index {
	my ($DBfna) = @_;
	my $lambdaIdxBin = configured_lambda3();
	die "Cannot build Lambda index: lambda3 is not configured or executable\n" unless ($lambdaIdxBin ne "");
	my $BlastCores = 8; #just pick reasonable number
	print "###################################\nCompiling lambda database for $DBfna using $BlastCores cores\n";
	run_cmd($lambdaIdxBin, "mkindexn", "-t", $BlastCores, "-d", $DBfna);
	my $index = "$DBfna.lba";
	die "Lambda did not create expected index $index\n" unless (-s $index);
	my $pigz = command_exists("pigz");
	my @compressCmd;
	if ($pigz){
		@compressCmd = ($pigz, "-p", $BlastCores, $index);
	} else {
		my $gzip = command_exists("gzip") or die "gzip is required to compress $index\n";
		@compressCmd = ($gzip, $index);
	}
	my $compressed = "$index.gz";
	my $backup = "$compressed.installer-backup.$$";
	if (-e $compressed){
		rename($compressed, $backup) or die "Cannot preserve existing index $compressed: $!\n";
	}
	my $compressStatus = system(@compressCmd);
	if ($compressStatus != 0 || !-s $compressed){
		unlink($compressed) if (-e $compressed);
		rename($backup, $compressed) if (-e $backup);
		die "Could not compress Lambda index $index (status $compressStatus)\n";
	}
	if (-e $backup){
		unlink($backup) or warn "Could not remove old index backup $backup: $!\n";
	}
	print "Compiled index for $DBfna\n\n";
}

sub getbeetax($){
	my ($aref) = @_;
	my @txt = @{$aref};
	print "Downloading bee specific database and taxonomy.\n";
	ensure_dir("$ddir/beeTax/");
	my $DB = "$ddir/beeTax/beeTax.fasta"; my $DBtax = "$ddir/beeTax/beeTax.txt";
	#getS2("http://5.196.17.195/pr2/download/representative_sequence_of_each_cluster/gb203_pr2_all_10_28_99p.fasta.tar.gz",$DB.".tar.gz");
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.fna",$DB);
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.txt",$DBtax);
	getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/beeTax.fasta.lba.gz","$DB.lba.gz") if ($downloadLmbdIdx);
	#parse_PR2($DB,$DBtax); #unlink ($DBtax.".pre");
	@txt = addInfoLtS("TAX_REFDB_BEE",$DB,\@txt,1);
	@txt = addInfoLtS("TAX_RANK_BEE",$DBtax,\@txt,1);
	buildIndex($DB);
	return (@txt);
}

sub getPR2db($){
	my ($aref) = @_;
	my @txt = @{$aref};
	print "Downloading PR2 99% clustered database.\n";
	ensure_dir("$ddir/PR2/");
	#my $DB = "$ddir/PR2/PR2_pack"; 
	my $DBtax = "$ddir/PR2_5.0_tax.txt";
	#getS2("http://5.196.17.195/pr2/download/representative_sequence_of_each_cluster/gb203_pr2_all_10_28_99p.fasta.tar.gz",$DB.".tar.gz");
	#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/gb203PR2.tar.gz",$DB.".tar.gz");
#	system "tar -xzf $DB.tar.gz -C $ddir/PR2;rm $DB.tar.gz";
	#getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/gb203_pr2_all_10_28_99p.fasta.lba.gz","$ddir/PR2/gb203_pr2_all_10_28_99p.fasta.lba.gz") if ($downloadLmbdIdx);

#https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.tax.gz

	#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/PR2/pr2_version_5.0.0_SSU_mothur.tax.gz",$DBtax.".gz");
	
	#my $DB = "$ddir/PR2_5.0_pre.fasta";
	my $DB = "$ddir/PR2_5.0.fasta";
	#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/PR2//pr2_version_5.0.0_SSU_mothur.fasta.gz",$DB.".gz");
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/PR2//pr2_version_5.0.0_SSU_UTAX.fasta.gz",$DB.".gz");
	gunzip_file("$DB.gz", $DB);
	parse_PR2($DB,$DBtax);
	#die "$DB,$DBtax\n";
	#parse_PR2($DB,$DBtax); #unlink ($DBtax.".pre");
	@txt = addInfoLtS("TAX_REFDB_PR2",$DB,\@txt,1);
	@txt = addInfoLtS("TAX_RANK_PR2",$DBtax,\@txt,1);
	
	
	
	buildIndex($DB);
	return (@txt);
}
sub getHITdb($){
	my ($aref) = @_;
	my @txt = @{$aref};
	print "Downloading HITdb April 2015 release..\n";
	ensure_dir("$ddir/HITdb/");
	my $DB = "$ddir/HITdb/HITdb_sequences.fasta"; my $DBtax = "$ddir/HITdb/HITdb_taxonomy.txt";
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_sequences.fna",$DB);
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_taxonomy_qiime.txt",$DBtax.".pre");
	getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/HITdb_sequences.fasta.lba.gz","$DB.lba.gz") if ($downloadLmbdIdx);
	parse_hitdb($DBtax.".pre",$DBtax); unlink ($DBtax.".pre");
	@txt = addInfoLtS("TAX_REFDB_HITdb",$DB,\@txt,1);
	@txt = addInfoLtS("TAX_RANK_HITdb",$DBtax,\@txt,1);
	buildIndex($DB);
	return (@txt);
}


sub getGG2($){
	my ($aref) = @_;
	my @txt = @{$aref};
	#greengenes ------------------------
	my $gg1 = "https://ksgp.earlham.ac.uk/downloads/greengenes2/GG2.2022.10.fasta.gz";
	my $gg2 = "https://ksgp.earlham.ac.uk/downloads/greengenes2/GG2.2022.10.tax.gz";
	my $DB = "$ddir/GG2.2022.10.fasta";
	#system("wget -O $DB.gz $gg1");
	print "Downloading GreenGenes2 2022 release..\n";
	getS2($gg1,"$DB.gz");
	getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/GG2.2022.10.fasta.lba.gz","$DB.lba.gz") if ($downloadLmbdIdx);
	sleep(10);
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("TAX_REFDB_GG2",$DB,\@txt,1);
	buildIndex($DB);
	$DB = "$ddir/GG2.2022.10.tax";
	#system("wget -O $DB.gz $gg2");
	getS2($gg2,"$DB.gz");
	sleep(3);
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("TAX_RANK_GG2",$DB,\@txt,1);
	return @txt;
}



sub getKSGP($){
	my ($aref) = @_;
	my @txt = @{$aref};
	

	my $DB = "$ddir/KSGPv4.0";
	print "Downloading KSGP v4.0 Jul 2026 release..\n";
	my $tarUTN = "$ddir/KSGPv4.0.gz";	my $tarUTNtax = "$ddir/KSGPv4.0.tax.gz";
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_v4.0.fasta.gz",$tarUTN);
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_plus2.tax.gz",$tarUTNtax);
	gunzip_file($tarUTN, "$DB.fasta"); gunzip_file($tarUTNtax, "$DB.tax");
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGPv4.0.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);



#	my $DB = "$ddir/KSGP_v3.1";unlink glob("${DB}*");
#	print "Downloading KSGP v3.1 2025 release..\n";
#	my $tarUTN = "$ddir/KSGPv3.1.gz";	my $tarUTNtax = "$ddir/KSGPv3.1.tax.gz";
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1/KSGP.fasta.gz",$tarUTN);
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1/KSGP.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1//KSGP_v3.1.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);

#	my $DB = "$ddir/KSGP_v2.0";unlink glob("${DB}*"); print "Downloading KSGP v3 2025 release..\n";
#	my $tarUTN = "$ddir/KSGPv3.gz";	my $tarUTNtax = "$ddir/KSGPv3.tax.gz";
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3/KSGP_v3.fasta.gz",$tarUTN);
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3/KSGP_v3.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");


#	print "Downloading KSGP 2024 release..\n";
#	my $DB = "$ddir/KSGP_v2.0";unlink glob("${DB}*"); my $tarUTN = "$ddir/KSGPv2.gz";	my $tarUTNtax = "$ddir/KSGPv2.tax.gz";
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2/KSGP_v2.fasta.gz",$tarUTN);
#	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2/KSGP_LCA_v2.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");
	
	
	#getS2("https://ksgp.earlham.ac.uk/downloads/v1.0/KSGP_v1.0.tar.gz",$tarUTN);
	#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2//KSGP_v2.0.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);
	#system "tar -xzf $tarUTN -C $ddir;rm -f $tarUTN";
	@txt = addInfoLtS("TAX_RANK_KSGP","$DB.tax",\@txt,1);
	@txt = addInfoLtS("TAX_REFDB_KSGP","$DB.fasta",\@txt,1);
	buildIndex("$DB.fasta");
	print "Added $DB.fasta and $DB.tax to lotus config.\n";
	return @txt;
}

sub getSLV($){
	my ($aref) = @_;
	my @txt = @{$aref};
	my $locSLBdl = 0;
	#SILVA -----------------------------------
	#TAX_REFDB_SLV  TAX_REFDB_SLV
	#changed to ver 119
	#changed to 123
	#changed to 128
	#changed to 132
#	my $baseSP = "http://www.arb-silva.de/fileadmin/silva_databases/release_123_1/Exports";
	my $SLVver = "138.1";
	#my $baseSP = "http://www.arb-silva.de/fileadmin/silva_databases/release_$SLVver/Exports";
	my $baseSP = "https://ftp.arb-silva.de/release_$SLVver/Exports";
#	my $baseSN = "SILVA_123.1";my $baseLN = "SLV_123.1";	my $SLVver = "123.1";
	my $baseSN = "SILVA_$SLVver";	my $baseLN = "SLV_$SLVver";	
	
	my $DB2 = "$ddir/$baseLN"."_SSU.tax";
	my $DB = "$ddir/$baseLN"."_SSU.fasta";
	print "Downloading SILVA SSU release $SLVver..\n";
	if ($locSLBdl){ #in case silva server doesn't work again..
		$baseSP = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/";
		#my $SlvAltFna = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_SSU.fasta.gz";
		#getS2($SlvAltFna,"$DB.gz");
		#system("gunzip -c $DB.gz > $DB;rm -f $DB.gz"); 
		#my $SlvAltTax = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_SSU.tax.gz";
		#getS2($SlvAltTax,"$DB2.gz");
		#system("gunzip -c $DB2.gz > $DB2;rm -f $DB2.gz"); 
	} 
	
	my $SLV = $baseSP."/".$baseSN."_SSURef_NR99_tax_silva.fasta.gz";
	getS2($SLV,"$DB.gz");
	getS2("https://ksgp.earlham.ac.uk/lambdaDBs/v3.0/SLV_138.1_SSU.fasta.lba.gz","$DB.lba.gz") if ($downloadLmbdIdx);
	#print "$SLV\n";
	gunzip_file("$DB.gz", "$ddir/SSUsilva.fasta"); 
	getS2($baseSP."/taxonomy/tax_slv_ssu_$SLVver.txt.gz","$ddir/SLVtaxSSU.csv.gz");
	#print "$baseSP/taxonomy/tax_slv_ssu_$SLVver.txt.gz\n";
	gunzip_file("$ddir/SLVtaxSSU.csv.gz", "$ddir/SLVtaxSSU.csv"); 
	prepareSILVA("$ddir/SSUsilva.fasta",$DB,$DB2,"$ddir/SLVtaxSSU.csv","");
	unlink("$ddir/SSUsilva.fasta");

	$finalWarning .= "\nWARNING: Silva $SLVver does not have consistent taxonomy levels for LSU's, therefore the taxonomy used in LotuS will contain \"?\" after taxonomy name.\n";
	
	@txt = addInfoLtS("TAX_REFDB_SSU_SLV",$DB,\@txt,1);
	@txt = addInfoLtS("TAX_RANK_SSU_SLV",$DB2,\@txt,1);
	buildIndex($DB);
	
#------------------------------ LSU SLV DB --------------------------
	$DB = "$ddir/$baseLN"."_LSU.fasta";
	$DB2 = "$ddir/$baseLN"."_LSU.tax";
	print "Downloading SILVA LSU release $SLVver..\n";
	$locSLBdl=0; $SLVver="138.1";#change this to local (132 release), since SIVLA doesn't have that yet..
	if ($locSLBdl){ #in case silva server doesn't work again..
		$baseSP = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/";
	}
	#	my $SlvAltFna = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_LSU.fasta.gz";
	#	getS2($SlvAltFna,"$DB.gz");
	#	system("gunzip -c $DB.gz > $DB;rm -f $DB.gz"); 
	#	my $SlvAltTax = "https://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_LSU.tax.gz";
	#	getS2($SlvAltTax,"$DB2.gz");
	#	system("gunzip -c $DB2.gz > $DB2;rm -f $DB2.gz"); 
	$SLV = $baseSP."/".$baseSN."_LSURef_tax_silva.fasta.gz";
	getS2($SLV,"$DB.gz");
	getS2($baseSP."/taxonomy/tax_slv_lsu_$SLVver.txt.gz","$ddir/SLVtaxLSU.csv.gz");
	gunzip_file("$ddir/SLVtaxLSU.csv.gz", "$ddir/SLVtaxLSU.csv");
	gunzip_file("$DB.gz", "$ddir/LSUSILVA.fasta"); #unlink("$DB.tgz");
	prepareSILVA("$ddir/LSUSILVA.fasta",$DB,$DB2,"$ddir/SLVtaxLSU.csv","$ddir/SLVtaxSSU.csv");
	unlink("$ddir/LSUSILVA.fasta"); unlink("$ddir/SLVtaxLSU.csv");unlink("$ddir/SLVtaxSSU.csv");
	@txt = addInfoLtS("TAX_REFDB_LSU_SLV",$DB,\@txt,1);
	@txt = addInfoLtS("TAX_RANK_LSU_SLV",$DB2,\@txt,1);
	buildIndex($DB);

	return @txt;
}

sub prepareSILVA($ $ $ $ $){
	#taxf3 is for 18S/28S #taxf3 is for SSU/LSU
	my ($path, $SeqF,$taxF,$taxGuide,$taxGuide2) = @_;
	print("Rewriting SILVA DB..\n");
	my %taxG;


	open I,"<",$taxGuide or die "Can't find taxguide file $taxGuide\n";
	while (my $line = <I>){
		chomp($line); my @splg = split("\t",$line);
		if (scalar(@splg) > 2){
			my $newN =  $splg[0]; #lc
			$taxG{$newN} =  $splg[2];
		}
	} 
	close I;

	if ($taxGuide2 ne ""){
	open I,"<",$taxGuide2 or die "Can't find taxguide file $taxGuide2\n";
	while (my $line = <I>){
		chomp($line); my @splg = split("\t",$line);
		next if (@splg < 3);
		my $newN =  $splg[0]; #lc
		$taxG{$newN} =  $splg[2];
	} 
	close I;
	}

	my $taxTmp = "$taxF.tmp.$$";
	my $seqTmp = "$SeqF.tmp.$$";
	open I,"<",$path or die ("could not find SILVA file \n$path\n");
	open OT,">",$taxTmp or die "Cannot write SILVA taxonomy $taxTmp: $!\n";
	open OS,">",$seqTmp or die "Cannot write SILVA sequences $seqTmp: $!\n";
	#open OT2,">",$taxF2;open OS2,">",$SeqF2;
	my @tdesign = (" k__"," p__"," c__"," o__"," f__"," g__"," s__");
	my $skip = 0;
	my $eukMode = 0;
	my $replacementTax =0; my $allTax=0;
	while (my $line = <I>){
		chomp($line);
		if ($line =~ m/^>/){#header
			$skip=0;$eukMode = 0;
			my @spl = split("\\.",$line);
			if (1){
				; #do nothing
			}elsif ($spl[0] =~ m/>AB201750/){
				$line = ">AB201750.1.1495 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 2;Anaerovirgula;Anaerovirgula multivorans";
				@spl = split("\\.",$line);
			} elsif ($spl[0] =~ m/>DQ643978/){
				$line = ">DQ643978.1.1627 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Geosporobacter;Geosporobacter subterraneus";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>X99238/){
				$line = ">X99238.1.1404 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 1;Thermobrachium;Thermobrachium celere";
				@spl = split("\\.",$line);
			} elsif ($spl[0] =~ m/>FJ481102/){
				$line = ">FJ481102.1.1423 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 1;Fervidicella;Fervidicella metallireducens AeB";
				@spl = split("\\.",$line);
			} elsif ($spl[0] =~ m/>EU443727/){
				$line = ">EU443727.1.1627 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;Thermotalea metallivorans";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>FR690973/){
				$line = ">FR690973.1.2373 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Thiopilula aggregata";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>CP002161/){
				$line = ">CP002161.5310.6845 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Candidatus Zinderia;Candidatus Zinderia insecticola CARI";
				@spl = split("\\.",$line);
			} elsif ($spl[0] =~ m/>FR690975/){
				$line = ">FR690975.1.2297 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Thiopilula aggregata";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>FR690991/){
				$line = ">FR690991.1.2147 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Marithioploca araucae";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>FR690991/){
				$line = ">AB910318.1.1553 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>AB910318/){
				$line = ">AB910318.1.1553 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
				@spl = split("\\.",$line);
			}elsif ($spl[0] =~ m/>AY796047/){
				$line = ">AY796047.1.1592 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
				@spl = split("\\.",$line);
			}

			
			my $ID = $spl[0];
			$ID = substr($ID,1);
			$line =~ m/[^\s]+\s(.*)$/;
			my $tax = $1; chomp $tax;
			if ($tax =~ m/^\s*Eukaryota/){$eukMode = 1;}#$skip = 1; next;}
			
			print OS ">".$ID."\n";
			@spl = split(";",$tax);
			for (my $i=0;$i<@spl;$i++){
				$spl[$i] =~ s/^\s*//; $spl[$i] =~ s/\s*$//;
			}
			#die "@spl\n";
			my $tline;
			if (!$eukMode){
				if (@spl > 7 ){
					print $line."\n";
					print("too many categories\n");
				}
				for (my $i=0;$i<7; $i++){
					if ($i < scalar(@spl)){
						if ($spl[$i] =~ m/^unidentified/){$spl[$i] = "?";}
						$spl[$i] = $tdesign[$i].$spl[$i];
					} else {
						$spl[$i] = $tdesign[$i];
					}
				}
				$tline = $ID ."\t".join(";",@spl);
			} else {#parse the levels out from taxguide
				my $tmpTax = "";
				my @jnd;
				my @soughtCls = ("domain","phylum","class","order","family","genus","species");
				my $soughtLvl = 0;  my $lastUsed = 0;
				for (my $i=0;$i<@spl; $i++){
					my $scanTax = $tmpTax.$spl[$i].";";
					if (exists($taxG{$scanTax}) || $soughtLvl == 6 || $spl[$i] =~ m/^unidentified/){
						#print "$taxG{$scanTax} LL\n";
						#SILVA has no species level in tax guide file
						$lastUsed = $soughtLvl;
						if ($soughtLvl == 6){
							push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
							$soughtLvl++;
							last;
						} elsif ($spl[$i] =~ m/^unidentified/ || $taxG{$scanTax} eq ""){#Euk in LSU file have no annotation..
							$spl[$i] = "";
							push(@jnd,$tdesign[$soughtLvl]."?");
							$soughtLvl++;
						} elsif ($taxG{$scanTax} eq $soughtCls[$soughtLvl]){
							push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
							#print $tdesign[$soughtLvl].$spl[$i]."\n";
							$soughtLvl++;
						} elsif ($taxG{$scanTax} eq $soughtCls[$soughtLvl+1]){#fill in empty levels
							push(@jnd,$tdesign[$soughtLvl]);
							$soughtLvl++;
							push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
							#print "Skipped to level ".$tdesign[$soughtLvl].$spl[$i]."\n";
							$soughtLvl++;
						}
						
					} else { #more likely to be low level species
						my $arS = @spl;
						#species signatuer & last entry
						if ($spl[$i] =~ m/\S+\s\S+/ && $arS >= ($i)){
							my $ncnt=1;
							while ($soughtLvl<6){
								my $nIdx = $lastUsed+ $ncnt;
								if ($nIdx < ($arS-1) ){
									#just impute preceding levels
									push(@jnd,$tdesign[$soughtLvl]."?".$spl[ $nIdx ]);
								} else {
									push(@jnd,$tdesign[$soughtLvl]."?");
								}
								$soughtLvl++;$ncnt++;
							}
							$soughtLvl = 6;
							#almost certainly a species
							push(@jnd,$tdesign[$soughtLvl].$spl[$i]); 
							$soughtLvl++;
							$replacementTax++;
							#print $ID."\t".join(";",@jnd)."\n$lastUsed\n";
							last;
						
						} else {
							#print $scanTax." JJ\n";
						}
					}
					 #Eukaryota;Fungi;Ascomycota;Archaeorhizomycetes;Archaeorhizomycetales;Archaeorhizomycetales_incertae_sedis
					$tmpTax .= $spl[$i].";";
					$lastUsed = $i;
				}
				$allTax++;
				for (;$soughtLvl<7;$soughtLvl++){
					push(@jnd,$tdesign[$soughtLvl]);
				}
				$tline = $ID."\t".join(";",@jnd);
				#die $tax." CC " .$tline."\n";
			}
			print OT $tline."\n";
			#die($tline);
		} elsif ($skip == 0){ #work through sequence
			$line =~ s/\s//g;
			$line =~ s/U/T/g;
			$line =~ s/u/t/g;
			#die $line;
			print OS $line."\n";
		}
	}
	#print "$replacementTax out of $allTax could not be defined to clear taxonomic levels and were imputed (with mostly empty tax levels or a \"?\" before tax name\n";

	close I; close OT or die "Cannot close $taxTmp: $!\n"; close OS or die "Cannot close $seqTmp: $!\n"; #close OT2; close OS2;
	die "SILVA conversion produced empty taxonomy output $taxTmp\n" unless (-s $taxTmp);
	die "SILVA conversion produced empty sequence output $seqTmp\n" unless (-s $seqTmp);
	rename($taxTmp, $taxF) or die "Cannot replace $taxF with $taxTmp: $!\n";
	rename($seqTmp, $SeqF) or die "Cannot replace $SeqF with $seqTmp: $!\n";
}



# Pinned checksum for a URL, or undef (the caller then has to skip the download).
sub pinned_sha256 {
	my ($url) = @_;
	return $PINNED_SHA256{$url};
}

# Download $in to $out and verify its SHA-256 before $out is created. $expected overrides
# the pinned table (the ONT tools carry their own checksums). A copy of the same file in
# bin/installs/ with the right checksum is used instead of downloading it again.
sub getS2($ $;$){
	my ($in,$out,$expected) = @_;
	print "getS2:$in\n$out\n";
	die "Refusing non-https download URL: $in\n" unless ($in =~ m{^https://}i);
	$expected //= pinned_sha256($in);
	die "No pinned SHA-256 checksum for $in; refusing to install an unverified download.\n"
		unless (defined($expected) && $expected =~ /^[0-9a-f]{64}$/i);
	ensure_dir(dirname($out));
	my $tmp = "$out.tmp.$$";
	unlink($tmp) if (-e $tmp);
	(my $bundled = $in) =~ s{^.*/}{$bdir/installs/};
	if (-s $bundled && eval { verify_sha256($bundled, $expected) }) {
		print "Using bundled copy $bundled\n";
		copy($bundled, $tmp) or die "Can't copy $bundled to $tmp: $!\n";
	} elsif ($WGETpres){
		print "wget -O $tmp $in\n";
		run_cmd("wget", "-O", $tmp, $in);
	} elsif ($CURLpres){
		#-f: an HTTP error is a failure, not an error page saved as the file
		run_cmd("curl", "-fsSL", "--retry", "3", "-o", $tmp, $in);
	} elsif (http_tiny_https()){
		print "HTTP::Tiny $in\n";
		my $res = HTTP::Tiny->new->mirror($in, $tmp);
		die "Download failed for $in: HTTP $res->{status} $res->{reason}\n" unless ($res->{success});
	} else {
		die "Downloads need \"wget\" or \"curl\" (or the Perl modules IO::Socket::SSL and Net::SSLeay for HTTP::Tiny). Please install wget or curl.\n";
	}
	die "Download produced no file: $in -> $tmp\n" unless (-e $tmp);
	die "Downloaded file is empty: $in -> $tmp\n" unless (-s $tmp);
	eval { verify_sha256($tmp, $expected); 1 } or do {
		my $err = $@;
		unlink($tmp) if (-e $tmp);
		die "Download of $in failed verification; nothing was installed.\n$err";
	};
	rename($tmp, $out) or do {
		unlink($tmp) if (-e $tmp);
		die "Can't replace $out with downloaded file $tmp: $!\n";
	};
	return $out;
}


sub version_is_newer {
	my ($candidate,$current) = @_;
	return 0 unless (defined($candidate) && defined($current));
	return 0 unless ($candidate =~ m/^\d+(?:\.\d+)+$/ && $current =~ m/^\d+(?:\.\d+)+$/);
	my @candidateParts = split(/\./,$candidate);
	my @currentParts = split(/\./,$current);
	my $parts = @candidateParts > @currentParts ? scalar(@candidateParts) : scalar(@currentParts);
	for (my $i=0; $i<$parts; $i++){
		my $candidatePart = $candidateParts[$i] // 0;
		my $currentPart = $currentParts[$i] // 0;
		return 1 if ($candidatePart > $currentPart);
		return 0 if ($candidatePart < $currentPart);
	}
	return 0;
}

sub compile_LCA($){
	my ($ldi2) = @_;
	my $expPath = "$bdir/LCA";
	if (-x $expPath){#test if can execute locally
		my ($lcaV,$status) = capture_cmd($expPath, "-v");
		return $expPath if ($status == 0 && $lcaV =~ m/0\.\d+/);
	}
	if (-d $ldi2 && -f "$ldi2/Makefile" ){
		print "Compiling LCA..\n";
		unlink bsd_glob("$ldi2/*.o");
		my $stat = system("make", "-C", $ldi2);
		if ($stat == 0){
			unlink("$ldir/LCA") if -e "$ldir/LCA"; unlink("$bdir/LCA") if -e "$bdir/LCA"; move("$ldi2/LCA", "$bdir/LCA") or die "Cannot install LCA: $!\n"; run_cmd("chmod", "+x", "$bdir/LCA");
		} else {
			die "Compilation of required LCA binary failed (make status $stat). Install a C++ compiler and rerun the installer.\n";
		}
	} else {
		die "LCA source directory or Makefile is missing at $ldi2\n";
	}
	die "Compilation did not produce executable $expPath\n" unless (-e $expPath);
	run_cmd("chmod", "+x", $expPath);
	return $expPath;
}
sub install_bundled_rtk($){
	my ($source) = @_;
	die "Bundled rtk binary is missing or empty at $source\n" unless (-s $source);
	run_cmd("chmod", "+x", $source);
	my ($rtk_help,$source_status) = capture_cmd($source, "-h");
	die "Bundled rtk binary at $source could not be executed or did not identify itself as rtk (help exit status $source_status).\n"
		unless ($rtk_help =~ /rarefaction tool kit \(rtk\)\s+[\d.]+/i);

	my $destination = "$bdir/rtk";
	my $source_abs = abs_path($source) // $source;
	my $destination_abs = -e $destination ? (abs_path($destination) // $destination) : $destination;
	return $destination if $source_abs eq $destination_abs;
	copy_file_atomic($source, $destination);
	run_cmd("chmod", "+x", $destination);
	my ($installed_help,$installed_status) = capture_cmd($destination, "-h");
	die "Installed rtk binary at $destination failed its execution check (help exit status $installed_status).\n"
		unless ($installed_help =~ /rarefaction tool kit \(rtk\)\s+[\d.]+/i);
	return $destination;
}

sub compile_sdm($){
	my ($ldi2) = @_;
	my $expPath = "$bdir/sdm";
	if (-x $expPath){#test if can execute locally
		my ($sdmV,$status) = capture_cmd($expPath, "-v");
		return $expPath if ($status == 0 && $sdmV =~ m/sdm \d/);
	}
	if (-d $ldi2 && -f "$ldi2/Makefile" && -f "$ldi2/DNAconsts.cpp"){
		print "Compiling sdm..\n";
		unlink bsd_glob("$ldi2/*.o");
		my $stat = system("make", "-C", $ldi2);
		if ($stat != 0){#repeat without gzip
			print "\n\n\n\n=================\nProblem compiling sdm with gzip support\nFallback to sdm compilation without gzip support\n";
			my $header = "$ldi2/DNAconsts.h";
			my $backup = "$header.installer-backup.$$";
			copy($header, $backup) or die "Cannot back up $header before fallback compilation: $!\n";
			run_cmd($^X, "-pi", "-e", "s/#define _gzipread/#define _notgzip/g", $header);
			unlink bsd_glob("$ldi2/*.o");
			$stat = system("make", "-C", $ldi2);
			copy($backup, $header) or die "Cannot restore $header after fallback compilation: $!\n";
			unlink($backup) or warn "Could not remove temporary backup $backup: $!\n";
			$finalWarning .= "Can not read gzip file\n";
		}
		if ($stat == 0){
			unlink("$ldir/sdm") if -e "$ldir/sdm"; unlink("$bdir/sdm") if -e "$bdir/sdm"; move("$ldi2/sdm", "$bdir/sdm") or die "Cannot install sdm: $!\n"; run_cmd("chmod", "+x", "$bdir/sdm");
		} else {
			die "Compilation of required sdm binary failed (make status $stat). Install a C++ compiler and rerun the installer.\n";
		}
	} else {
		die "sdm source directory or required source files are missing at $ldi2\n";
	}
	die "Compilation did not produce executable $expPath\n" unless (-e $expPath);
	run_cmd("chmod", "+x", $expPath);
	return $expPath;
}

sub command_exists {
	my ($cmd) = @_;
	return unless defined($cmd) && $cmd =~ m/^[A-Za-z0-9_.+\-]+$/;
	foreach my $dir (split(/:/, $ENV{PATH} // "")){
		my $path = "$dir/$cmd";
		return $path if (-x $path);
	}
	return;
}

sub run_cmd {
	my (@cmd) = @_;
	die "run_cmd called without command\n" unless @cmd;
	print "+ @cmd\n";
	my $status = system(@cmd);
	die "Could not execute $cmd[0]: $!\n" if ($status == -1);
	die "Command failed: @cmd\nExit status: ".($status >> 8)."\n" if ($status != 0);
}

sub capture_cmd {
	my (@cmd) = @_;
	die "capture_cmd called without command\n" unless @cmd;
	print "+ @cmd\n";
	open(my $fh, "-|", @cmd) or die "Could not execute $cmd[0]: $!\n";
	local $/;
	my $output = <$fh>;
	$output = "" unless defined($output);
	close($fh);
	my $raw_status = $?;
	my $status = $raw_status == -1 ? -1 : ($raw_status & 127) ? 128 + ($raw_status & 127) : ($raw_status >> 8);
	return ($output,$status);
}

sub read_user_input {
	my ($context, $choices, $default) = @_;
	while (1) {
		my $line = <STDIN>;
		die "End of input while waiting for $context; installation aborted.\n" unless defined($line);
		chomp($line);
		return $line unless defined($choices);
		$line =~ s/^\s+|\s+$//g;
		$line = $default if $line eq "" && defined($default);
		return $line if grep { $line eq $_ } @$choices;
		print "Invalid answer; enter " . join(" or ", @$choices) . ": ";
	}
}

sub ensure_dir {
	my ($dir) = @_;
	die "Refusing empty directory path\n" unless defined($dir) && $dir ne "";
	make_path($dir) unless (-d $dir);
}

sub gunzip_file {
	my ($in,$out) = @_;
	die "Missing gzip input $in\n" unless (-s $in);
	my $tmp = "$out.tmp.$$";
	unlink($tmp) if (-e $tmp);
	if (!gunzip $in => $tmp){
		unlink($tmp) if (-e $tmp);
		die "gunzip failed for $in -> $tmp: $GunzipError\n";
	}
	if (!-s $tmp){
		unlink($tmp) if (-e $tmp);
		die "gunzip produced empty output $tmp\n";
	}
	rename($tmp, $out) or do {
		unlink($tmp) if (-e $tmp);
		die "Cannot replace $out with decompressed file $tmp: $!\n";
	};
	unlink($in) or warn "Could not remove $in: $!\n";
}

sub write_config_atomic {
	my ($cfg,$lines) = @_;
	my $tmp = "$cfg.tmp.$$";
	if (-e $cfg && !$configBackupWritten){
		copy_file_atomic($cfg, "$cfg.bak");
		$configBackupWritten = 1;
	}
	open(my $fh, ">", $tmp) or die "Cannot write $tmp: $!\n";
	print {$fh} @{$lines} or die "Cannot write $tmp: $!\n";
	close($fh) or die "Cannot close $tmp: $!\n";
	rename($tmp, $cfg) or die "Cannot replace $cfg with $tmp: $!\n";
}

sub copy_file_atomic {
	my ($source,$destination) = @_;
	die "Cannot copy missing or empty file $source\n" unless (-s $source);
	my $tmp = "$destination.tmp.$$";
	unlink($tmp) if (-e $tmp);
	copy($source, $tmp) or die "Cannot copy $source to $tmp: $!\n";
	die "Copy of $source to $tmp is incomplete\n" unless (-s $tmp == -s $source);
	rename($tmp, $destination) or do {
		unlink($tmp) if (-e $tmp);
		die "Cannot replace $destination with $tmp: $!\n";
	};
}

sub replace_tree_atomic {
	my ($source,$destination) = @_;
	die "Cannot install missing directory $source\n" unless (-d $source);
	my $backup = "$destination.installer-backup.$$";
	die "Refusing to overwrite stale update backup $backup\n" if (-e $backup);
	if (-d $destination){
		rename($destination, $backup) or die "Cannot back up $destination to $backup: $!\n";
	}
	if (!rename($source, $destination)){
		my $error = $!;
		rename($backup, $destination) if (-d $backup);
		die "Cannot install $source as $destination: $error\n";
	}
	remove_tree($backup) if (-d $backup);
}

sub verify_sha256 {
	my ($file,$expected) = @_;
	return 1 unless defined($expected) && $expected ne "";
	require Digest::SHA;
	open(my $fh, "<:raw", $file) or die "Cannot open $file for checksum: $!\n";
	my $got = Digest::SHA->new(256)->addfile($fh)->hexdigest;
	close($fh);
	die "Checksum mismatch for $file\nExpected: $expected\nGot: $got\n" unless (lc($got) eq lc($expected));
	return 1;
}

sub check_version {
	my ($cmd) = @_;
	my $exe = (-x $cmd) ? $cmd : command_exists($cmd);
	return 0 unless $exe;
	#R before 4.2 prints "Rscript --version" to stderr
	my ($check,$status) = capture_cmd_merged($exe, "--version");
	return 0 if ($status != 0 && $check eq "");
	if ($check =~ m/version\s+([0-9]+)(?:\.[0-9]+)*/){
		return $1;
	}
	if ($check =~ m/\b([0-9]+)(?:\.[0-9]+)+\b/){
		return $1;
	}
	return 0;
}

sub getTaxSfromUNITE{
	my ($head) = @_;
	my $taxS = "k__?;p__?;c__?;o__?;f__?;g__?;s__";
	if ($head =~ s/\|([^\|]+)$//){
		$taxS = $1;
	} else {
		die "Error in extrTaxFromFasta:: can't find \"|\" in string $head\n";
	}
	return ($taxS, $head);
}

sub extrTaxFromFasta($ $ $){
	my ($inFA, $oFA, $oTax) = @_;
	my $fastaTmp = "$oFA.tmp.$$";
	my $taxTmp = "$oTax.tmp.$$";
	open I,"<$inFA" or die "Can't open inFA $inFA\n";
	open OF,">",$fastaTmp or die "Can't write FASTA output $fastaTmp: $!\n";
	open OT,">",$taxTmp or die "Can't write taxonomy output $taxTmp: $!\n";
	
	my $fasta="";my $head="";#my $line="";
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^>/ ){
			if ($head ne ""){
				my ($taxS,$h2) = getTaxSfromUNITE($head);
				#die "$taxS\n$head\n";
				print OF ">$h2\n$fasta\n";
				print OT "$h2\t$taxS\n";
			}
			$fasta = "";
			$head = substr($line,1);
			next;
		}
		$fasta .= $line;
	}
	#final round..
	my ($taxS,$h2) =getTaxSfromUNITE($head);
	print OF ">$h2\n$fasta\n";
	print OT "$h2\t$taxS\n";
	
	
	close I; close OF or die "Cannot close $fastaTmp: $!\n"; close OT or die "Cannot close $taxTmp: $!\n";
	die "UNITE conversion produced empty FASTA output $fastaTmp\n" unless (-s $fastaTmp);
	die "UNITE conversion produced empty taxonomy output $taxTmp\n" unless (-s $taxTmp);
	rename($fastaTmp, $oFA) or die "Cannot replace $oFA with $fastaTmp: $!\n";
	rename($taxTmp, $oTax) or die "Cannot replace $oTax with $taxTmp: $!\n";
}


sub get_DBs{
#-------BIG DB INSTALL
	if ($refDBinstall[2] || $refDBinstall[8]){
		@txt = getSLV(\@txt);
	}
	if ($refDBinstall[1] || $refDBinstall[8]){
		@txt = getKSGP(\@txt);
	}
	if ($refDBinstall[3] || $refDBinstall [8]){
		@txt = getGG2(\@txt);
	}

	if ($refDBinstall [4] || $refDBinstall [8]){
		@txt = getHITdb(\@txt);
	}
	if ($refDBinstall [5] || $refDBinstall [8]){
		@txt = getPR2db(\@txt);
	}
	if ($refDBinstall [6] || $refDBinstall [8]){
		@txt = getbeetax(\@txt);
	}

	if ($refDBinstall[0]){
		print "No Ref DB will be installed.\n";
	}

	if ($getUTAX){
		print "Downloading UTAX ref databases..\n";
		my $tarUTN = "$ddir/utax_16s.tar.gz";
		getS2("https://drive5.com/utax/data/utax_rdp_16s_tainset15.tar.gz",$tarUTN);
		run_cmd("tar", "-xzf", $tarUTN, "-C", $ddir); unlink($tarUTN) or warn "Could not remove $tarUTN: $!\n";
		$tarUTN="$ddir/utax_ITS.tar.gz";
		getS2("https://drive5.com/utax/data/utax_unite_v7.tar.gz",$tarUTN);
		run_cmd("tar", "-xzf", $tarUTN, "-C", $ddir); unlink($tarUTN) or warn "Could not remove $tarUTN: $!\n";
		@txt = addInfoLtS("TAX_REFDB_SSU_UTAX","$ddir/utaxref/rdp_16s_trainset15/",\@txt,2);
		@txt = addInfoLtS("TAX_REFDB_ITS_UTAX","$ddir/utaxref/unite_v7/",\@txt,2);
		#die "X\n";
		
	}
	
	

	
	#-------BIG DB INSTALL END
	
	
	if ($ITSready){
		#ITS DB
		#my $tarUN = "$ddir/qITSfa.zip";
		#v9 2023 releast
		#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_refs_qiime_ver8_99_s_all_02.02.2019.fasta.zip",$tarUN);
#		my $UNITEdb = "$ddir/UNITE/sh_refs_qiime_ver8_99_s_all_02.02.2019.fasta";
		#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/sh_qiime_release_02.03.2015.zip",$tarUN);
		#system("rm -fr $ddir/UNITE;unzip -q -o $tarUN -d $ddir/UNITE/");
		#getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_taxonomy_qiime_ver8_99_s_all_02.02.2019.txt.zip",$tarUN);
		#system("unzip -q -o $tarUN -d $ddir/UNITE/;rm -rf $ddir/UNITE/__MACOSX/");
		#@txt = addInfoLtS("TAX_RANK_ITS_UNITE","$ddir/UNITE/sh_taxonomy_qiime_ver8_99_s_all_02.02.2019.txt",\@txt,1);
		my $tarUN = "$ddir/qITSfa.gz";
		#my $dlUNITE = "https://lotus2.earlham.ac.uk/lotus/packs/UNITE/v9_Dec23/sh_general_release_dynamic_all_25.07.2023.fasta.gz";
		my $dlUNITE = "https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_general_release_dynamic_s_all_19.02.2025.fasta.gz";
		getS2($dlUNITE,$tarUN);

#		my $UNITEdb = "$ddir/UNITE/sh_refs_v9_25.07.2023";
		my $UNITEdb = "$ddir/UNITE/sh_refs_v10_19.02.2025";
		
		ensure_dir("$ddir/UNITE/"); gunzip_file($tarUN, "$UNITEdb.fasta.tmp");
		extrTaxFromFasta("$UNITEdb.fasta.tmp","$UNITEdb.fasta","$UNITEdb.tax");
		unlink("$UNITEdb.fasta.tmp") or warn "Could not remove $UNITEdb.fasta.tmp: $!\n";
		
		
		#index creation/download
		#getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/sh_refs_v9_25.07.2023.fasta.lba.gz","$UNITEdb.fasta.lba.gz") if ($downloadLmbdIdx);
		getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/Lambda3/sh_refs_v10_19.02.2025.fasta.lba.gz","$UNITEdb.fasta.lba.gz") if ($downloadLmbdIdx);
		buildIndex("$UNITEdb.fasta");
		
		@txt = addInfoLtS("TAX_REFDB_ITS_UNITE","$UNITEdb.fasta",\@txt,1);
		@txt = addInfoLtS("TAX_RANK_ITS_UNITE","$UNITEdb.tax",\@txt,1);

		unlink($tarUN);
	}
	
	if ($ITSready){#ITS chimera check ref DB
		#my $itsDB = "https://lotus2.earlham.ac.uk/lotus/packs/DB/uchime_reference_dataset_11.03.2015.zip";
		my $itsDB = "https://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/uchime/uchime_UNITE_16_10_22.zip";
		getS2($itsDB,"$ddir/uchITS.zip");
		my $uchimeD = "$ddir/ITS_chimera22/";
		ensure_dir($uchimeD);
		run_cmd("unzip", "-q", "-o", "$ddir/uchITS.zip", "-d", $uchimeD);
		unlink("$ddir/uchITS.zip");
		$uchimeD .= "/2022_10_26_chimera_reference_release/";
		#die "$uchimeD/uchime_sh_refs_dynamic_original_985_11.03.2015.fasta";
		@txt = addInfoLtS("UCHIME_REFDB_ITS","$uchimeD/uchime_reference_dataset_16_10_2022.fasta",\@txt,1);
#		@txt = addInfoLtS("UCHIME_REFDB_ITS1","$uchimeD/ITS1_ITS2_datasets/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta",\@txt,1);
#		@txt = addInfoLtS("UCHIME_REFDB_ITS2","$uchimeD/ITS1_ITS2_datasets/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta",\@txt,1);
		@txt = addInfoLtS("UCHIME_REFDB_ITS1","$uchimeD/ITS1_ITS2_datasets/uchime_reference_dataset_16_20_2022_ITS1.fasta",\@txt,1);
		@txt = addInfoLtS("UCHIME_REFDB_ITS2","$uchimeD/ITS1_ITS2_datasets/uchime_reference_dataset_16_20_2022_ITS2.fasta",\@txt,1);
	}

	#-------------- install chimera check DBs
	
	# phiX ref genome
	my $phiXf = "$ddir/phiX.fasta";
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/DB/phiX.fasta",$phiXf);
	@txt = addInfoLtS("REFDB_PHIX",$phiXf,\@txt,1);
	
	#db gold #exchanged for rdp_gold since 1.30
	#my $goldDB = "http://drive5.com/uchime/gold.fa";
	my $goldDB = "https://lotus2.earlham.ac.uk/lotus/packs/rdp_gold.fa.gz";
	my $DB = "$ddir/rdp_gold.fa";
	#system("wget -O $DB $goldDB");
	getS2($goldDB,$DB.".gz");
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("UCHIME_REFDB",$DB,\@txt,1);



	#db Silva 119 clustered to 93% for LSUs
	my $LTUrefDB = "https://lotus2.earlham.ac.uk/lotus/packs/SILVA_119_LSU_93.ref.fasta.gz";
	$DB = "$ddir/SLV_119_LSU.fa";
	getS2($LTUrefDB,$DB.".gz");
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("UCHIME_REFDB_LSU",$DB,\@txt,1);
	
	
}
 

sub getUsearch{
	if ($isMac){
		my $usearch_path = command_exists("usearch");
		if ($usearch_path){
			print "Using USEARCH found at $usearch_path\n";
			@txt = addInfoLtS("usearch",$usearch_path,\@txt,1);
		} else {
			my $message = "No macOS USEARCH binary was found. Install USEARCH and rerun with -link_usearch /absolute/path/to/usearch.\n";
			print $message;
			$finalWarning .= $message;
		}
	} else {
		print "Downloading USEARCH v12 for sequence clustering and tax annotations..\n";
		my $usearch_path = "$bdir/usearch12_linux_beta";
		getS2("https://github.com/rcedgar/usearch12/releases/download/v12.0-beta1/usearch_linux_x86_12.0-beta",$usearch_path);
		run_cmd("chmod", "+x", $usearch_path);
		@txt = addInfoLtS("usearch",$usearch_path,\@txt,1);
	}
}


# Pin upstream releases so rerunning the installer does not silently change tools.
sub ont_tool_spec {
    my ($name) = @_;
    my %spec = (
        savont => ['0.7.0', 'bluenote-1577/savont', 'a60141fd4d4e83cdcbc3601220dc41487f4cdf7bd81558c7125e876117ed756c'],
        barbell => ['0.3.2', 'rickbeeloo/barbell', '2f840c3fe625c62f3d91deadf721da54a8c7b57f8f12a7976cc8cbba273f0daf'],
        minimap2 => ['2.28', 'lh3/minimap2', '5ea6683b4184b5c49f6dbaef2bc5b66155e405888a0790d1b21fd3c93e474278'],
    );
    die "Unknown ONT program $name\n" unless exists $spec{$name};
    return @{$spec{$name}};
}

sub ont_architecture {
    my @host = uname();
    my $arch = lc($host[4]);
    $arch = 'aarch64' if $arch eq 'arm64';
    $arch = 'x86_64' if $arch eq 'amd64';
    return $arch;
}

sub ont_rust_available {
    return 0 unless command_exists('cargo') && command_exists('rustc');
    my ($version, $status) = eval { capture_cmd(command_exists('rustc'), '--version') };
    return 0 if $@ || !defined($status) || $status != 0;
    my ($v) = $version =~ /rustc\s+(\d+\.\d+\.\d+)/;
    return defined($v) && !version_is_newer('1.88.0', $v);
}

sub savont_bioconda_binary {
    my $target = ($isMac ? 'osx-' : 'linux-') . ont_architecture();
    # Exact Bioconda 0.7.0 build-0 packages, including platform-specific digests.
    my %packages = (
        'linux-x86_64' => ['linux-64', 'hec9b1f2_0', 'e7ea28b084d176379d9fa273a3cb349c9e58e54436e2efd02c295acf91edbdd7'],
        'linux-aarch64' => ['linux-aarch64', 'h2013a2e_0', '6481120075ab330342dc15f856ae5a7cf2833d096275e53ec2af9189978bd80b'],
        'osx-x86_64' => ['osx-64', 'h121cdbd_0', '5270bb6e38f9eb315d88bd574df0bc0259d1cc3cb764a901f0197a53908feb01'],
        'osx-aarch64' => ['osx-arm64', 'ha819e4a_0', '3be62545d8ee9fa3b8ecd2e0b1c2aa4785075c915190468f71888ef571e3242f'],
    );
    die "No pinned Savont Bioconda binary is available for $target.\n" unless exists $packages{$target};
    my ($subdir, $build, $sha) = @{$packages{$target}};
    my ($version) = ont_tool_spec('savont');
    return ("https://api.anaconda.org/download/bioconda/savont/$version/$subdir/savont-$version-$build.conda", $sha);
}

sub savont_binary_failure {
    my ($reason) = @_;
    die "Savont Bioconda binary fallback failed:\n$reason\n"
        . "Install Rust >= 1.88 (including Cargo), a C/C++ compiler and CMake, then rerun the installer to compile Savont from source.\n"
        . "The existing savont configuration entry was preserved.\n";
}

sub ont_release_binary {
    my ($name) = @_;
    my $arch = ont_architecture();
    if ($name eq 'minimap2' && !$isMac && $arch eq 'x86_64') {
        return ('https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2',
            '51f2cf0e486d0f9f88ace1aa58fdc56571382a676ea0889ae607301c60693377', 'minimap2-2.28_x64-linux/minimap2');
    }
    if ($name eq 'barbell') {
        my $target = $arch . ($isMac ? '-apple-darwin' : '-unknown-linux-gnu');
        my %sha = (
            'aarch64-apple-darwin' => '2934846dcb2a6b2a6ebac86e67cf28b507b8b70a4440fd1f47213d3cf3181b87',
            'aarch64-unknown-linux-gnu' => '055a9fd1df729c671a8061ad26baa64fe078fe8a97e233b96812bea61977e9c1',
            'x86_64-unknown-linux-gnu' => 'f23a599eceb8b27211178facf5504935a98880c031eb9f2b4f136382193f2081',
        );
        return ("https://github.com/rickbeeloo/barbell/releases/download/v0.3.2/barbell-$target", $sha{$target}, '')
            if exists $sha{$target};
    }
    return; # Build from the pinned source archive on other supported architectures.
}

sub compatible_ont_program {
    my ($name, $path) = @_;
    return 0 unless defined($path) && -f $path && -x $path;
    my ($version, $status) = eval { capture_cmd($path, '--version') };
    return 0 if $@ || !defined($status) || $status != 0;
    if ($name eq 'minimap2') {
        my ($v) = $version =~ /(\d+\.\d+)/;
        return defined($v) && !version_is_newer('2.17', $v);
    }
    return 0 unless $version =~ /\b\Q$name\E\b/i;
    my ($help, $help_status) = eval { capture_cmd($path, $name eq 'savont' ? 'asv' : 'kit', '--help') };
    return 0 if $@ || !defined($help_status) || $help_status != 0;
    my @flags = $name eq 'savont'
        ? qw(--quality-value-cutoff --minimum-base-quality --chimera-allowable-errors --single-strand)
        : qw(--kit --input --output --maximize --threads);
    return 0 if grep { index($help, $_) < 0 } @flags;
    return 1;
}

sub find_ont_program {
    my ($name) = @_;
    my $configured = getInfoLtS($name, \@txt, '');
    $configured =~ s/"//g;
    my @candidates;
    if ($configured ne '' && $configured ne '??') {
        push @candidates, File::Spec->file_name_is_absolute($configured)
            ? $configured : File::Spec->catfile($ldir, $configured);
        push @candidates, $configured if -f $configured;
        push @candidates, command_exists($configured);
    }
    push @candidates, "$bdir/$name", command_exists($name);
    my %seen;
    for my $candidate (@candidates) {
        next unless defined($candidate) && $candidate ne '';
        my $path = abs_path($candidate);
        next unless defined($path) && !$seen{$path}++;
        return $path if compatible_ont_program($name, $path);
    }
    return;
}

sub ont_programs_to_install {
    # Minimap2 is also the default mapper for non-ONT workflows.
    return ('minimap2', ($installONT ? ('savont', 'barbell') : ()));
}

sub check_ont_build_requirements {
    my $needs_rust = 0; my $needs_make = 0;
    for my $name (ont_programs_to_install()) {
        next if find_ont_program($name);
        if ($name eq 'savont' && !ont_rust_available()) {
            eval {
                savont_bioconda_binary(); #fail early on unsupported platforms
                run_cmd($^X, "$ldir/helpers/extract_conda_executable.pl", '--check');
                1;
            } or savont_binary_failure($@);
            print "Usable Rust/Cargo not found; Savont will use the prebuilt Bioconda package.\n";
            next;
        }
        my @release = ont_release_binary($name);
        next if @release;
        if ($name eq 'minimap2') { $needs_make = 1; } else { $needs_rust = 1; }
    }
    if ($needs_rust) {
        for my $tool (qw(cargo rustc cc c++ cmake)) {
            die "ONT source builds require $tool. Install Rust >= 1.88, a C/C++ compiler and CMake, or install compatible Savont/Barbell executables on PATH, then rerun the installer.\n"
                unless command_exists($tool);
        }
        my ($version, $status) = capture_cmd(command_exists('rustc'), '--version');
        my ($rust_version) = $version =~ /rustc\s+(\d+\.\d+\.\d+)/;
        die "ONT source builds require Rust >= 1.88 (found: $version).\n"
            unless $status == 0 && defined($rust_version) && !version_is_newer('1.88.0', $rust_version);
    }
    if ($needs_make) {
        die "Building minimap2 requires make and a C compiler (plus zlib development headers).\n"
            unless command_exists('make') && command_exists('cc');
    }
}

# Everything a full install needs besides the ONT tools, checked before the database downloads.
sub check_full_install_requirements {
    my @missing;
    for my $tool (qw(tar gzip unzip make)) {
        push @missing, $tool unless command_exists($tool);
    }
    push @missing, 'a C compiler (gcc or cc)' unless command_exists('gcc') || command_exists('cc');
    push @missing, 'xz (unpacks the Lambda .tar.xz release)'
        if !$isMac && ($installBlast == 2 || $installBlast == 3) && !command_exists('xz');
    # sdm and LCA: the bundled Linux x86-64 executables, or their sources to compile
    for my $req (['sdm', qr/sdm \d/, 'Makefile'], ['LCA', qr/0\.\d+/, 'Makefile']) {
        my ($name, $versionRx, $makefile) = @$req;
        my $exe = "$bdir/$name";
        my $runs = 0;
        if (-f $exe) {
            chmod(0755, $exe);
            my ($version, $status) = capture_cmd_merged($exe, '-v');
            $runs = $status == 0 && $version =~ $versionRx;
        }
        push @missing, "a working $name: $exe does not run on this system and $ldir/${name}_src/$makefile is not present "
            . "(the bundled $name is a Linux x86-64 build; on other systems compile it from its source and place it at $exe)"
            unless $runs || -f "$ldir/${name}_src/$makefile";
    }
    die "The full installation needs:\n" . join("", map { " - $_\n" } @missing)
        . "Install these and rerun perl helpers/autoInstall.pl. Nothing has been downloaded yet.\n" if @missing;
    if (!command_exists('java')) {
        my $msg = "Java was not found: RDP taxonomy classification (-taxAligner 0) needs it at run time.\n";
        print $msg; $finalWarning .= $msg;
    }
}

sub install_ont_program {
    my ($name) = @_;
    if (my $existing = find_ont_program($name)) { return $existing; }
    my $bioconda = $name eq 'savont' && !ont_rust_available();
    my $installed = eval {
        my ($version, $repo, $sha) = ont_tool_spec($name);
        my $stage = tempdir("$name-install-XXXXXXXX", DIR => $bdir, CLEANUP => 1);
        my @release = ont_release_binary($name);
        my $exe;
        if ($bioconda) {
            my ($url, $digest) = savont_bioconda_binary();
            my $archive = "$stage/savont.conda";
            getS2($url, $archive, $digest);
            verify_sha256($archive, $digest);
            $exe = "$stage/savont";
            run_cmd($^X, "$ldir/helpers/extract_conda_executable.pl", $archive, 'bin/savont', $exe);
        } elsif (@release) {
            my ($url, $digest, $member) = @release;
            my $archive = "$stage/download";
            getS2($url, $archive, $digest);
            verify_sha256($archive, $digest);
            if ($member ne '') {
                run_cmd('tar', '-xjf', $archive, '-C', $stage);
                $exe = "$stage/$member";
            } else { $exe = $archive; }
        } else {
            my $archive = "$stage/source.tar.gz";
            getS2("https://codeload.github.com/$repo/tar.gz/refs/tags/v$version", $archive, $sha);
            verify_sha256($archive, $sha);
            run_cmd('tar', '-xzf', $archive, '-C', $stage);
            my $source = "$stage/$name-$version";
            if ($name eq 'minimap2') {
                my @host = uname();
                my @make = ('make', '-C', $source);
                push @make, 'arm_neon=1', 'aarch64=1' if $host[4] =~ /^(?:arm64|aarch64)$/;
                run_cmd(@make);
                $exe = "$source/minimap2";
            } else {
                my @build = (command_exists('cargo'), 'build', '--release', '--locked', '--manifest-path', "$source/Cargo.toml", '--target-dir', "$stage/target");
                push @build, '--config', "$source/.cargo/config.toml" if -f "$source/.cargo/config.toml";
                run_cmd(@build);
                $exe = "$stage/target/release/$name";
            }
        }
        die "Installation did not produce $name at $exe\n" unless -s $exe;
        run_cmd('chmod', '+x', $exe);
        die "Installed $name is not executable or lacks the CLI required by LotuS. See the tool's output above. The existing $name configuration entry was preserved.\n"
            unless compatible_ont_program($name, $exe);
        my $destination = "$bdir/$name";
        copy_file_atomic($exe, $destination);
        run_cmd('chmod', '+x', $destination);
        abs_path($destination);
    };
    if (my $error = $@) {
        savont_binary_failure($error) if $bioconda;
        die $error;
    }
    return $installed;
}

sub install_ont_programs {
    for my $name (ont_programs_to_install()) {
        my $path = install_ont_program($name);
        @txt = addInfoLtS($name, $path, \@txt, 1);
    }
}


# Run a command and return (stdout+stderr, exit status); some tools (vsearch, R < 4.2)
# print their version to stderr only.
sub capture_cmd_merged {
	my (@cmd) = @_;
	die "capture_cmd_merged called without command\n" unless @cmd;
	print "+ @cmd\n";
	my ($in, $out);
	my $pid = eval { open3($in, $out, undef, @cmd) }; #undef error handle: stderr joins stdout
	return ("", -1) unless $pid;
	close($in);
	my $output = do { local $/; <$out> } // "";
	close($out);
	waitpid($pid, 0);
	my $raw_status = $?;
	my $status = $raw_status == -1 ? -1 : ($raw_status & 127) ? 128 + ($raw_status & 127) : ($raw_status >> 8);
	return ($output, $status);
}

sub vsearch_works {
	my ($exe) = @_;
	return 0 unless defined($exe) && -f $exe;
	chmod(0755, $exe);
	my ($version, $status) = capture_cmd_merged($exe, "--version");
	return 0 unless $status == 0 && $version =~ m/vsearch v2\.(\d+)/;
	return $1 >= 15;
}

sub install_vsearch {
	my $bundled = "$bdir/vsearch";
	if (vsearch_works($bundled)) {
		print "Using the bundled vsearch at $bundled\n";
		@txt = addInfoLtS("vsearch", abs_path($bundled) // $bundled, \@txt, 1);
		return;
	}
	my $arch = ont_architecture();
	my $name = $isMac ? "vsearch-2.32.0-macos-universal"
		: $arch eq 'aarch64' ? "vsearch-2.32.0-linux-aarch64-static"
		: $arch eq 'x86_64' ? "vsearch-2.32.0-linux-x86_64" : "";
	if ($name eq "") {
		my $msg = "vsearch was not installed: no pinned vsearch release for architecture $arch (fallback to usearch).\n";
		print $msg; $finalWarning .= $msg;
		return;
	}
	print "Downloading vsearch 2.32.0 ($name)..\n";
	my $vtars = "$bdir/vsearch.tar.gz";
	getS2("https://github.com/torognes/vsearch/releases/download/v2.32.0/$name.tar.gz", $vtars);
	run_cmd("tar", "-xzf", $vtars, "-C", $bdir);
	unlink($vtars) or warn "Could not remove $vtars: $!\n";
	my $vexe = "$bdir/$name/bin/vsearch";
	if (vsearch_works($vexe)) {
		@txt = addInfoLtS("vsearch", $vexe, \@txt, 1);
	} else {
		my $msg = "vsearch at $vexe did not run on this system, so vsearch was not installed (fallback to usearch).\n";
		print "\n\nWARNING::\n$msg\n"; $finalWarning .= $msg;
	}
}

sub get_programs{
	my ($dtar, $dexe);
	#-----------  exit prog here, if set
	#-----------------------



	install_ont_programs();


	if ($ITSready){ #ITSx
		#itsx
		print "Downloading ITSX to detect valid ITS regions..\n";
		my $tarUTN = "$bdir/ITSx_1.1.4.tar.gz";
		getS2("https://lotus2.earlham.ac.uk/lotus/packs/ITSx_1.1.4.tar.gz",$tarUTN);
		run_cmd("tar", "-xzf", $tarUTN, "-C", $bdir); unlink($tarUTN) or warn "Could not remove $tarUTN: $!\n";
		@txt = addInfoLtS("itsx","$bdir/ITSx_1.1.4/./ITSx",\@txt,1);
		@txt = addInfoLtS("hmmsearch","$bdir/ITSx_1.1.4/bin/hmmsearch",\@txt,1);

	}

	#-------BLAST LAMBDA INSTALL
	if ($installBlast == 1 || $installBlast == 3){
		#Blast
		print "Downloading blast executables...\n";
		my $blfil = "ncbi-blast-2.2.29+-x64-linux.tar.gz";
		$exe = "$bdir/blast.tar.gz";
		if ($isMac){
			#the former macOS archive is no longer available upstream (HTTP 404)
			my $msg = "BLAST was not installed: no verified macOS BLAST+ package is available. Install blastn and makeblastdb (e.g. via conda or Homebrew) and set \"blastn\" and \"makeBlastDB\" in lOTUs.cfg, or use Lambda.\n";
			print $msg; $finalWarning .= $msg;
		} else {
		getS2("https://lotus2.earlham.ac.uk/lotus/packs/".$blfil,$exe);

		#my $path = "blast/executables/blast+/2.2.29/";
		#my $host = "ftp.ncbi.nlm.nih.gov";my $ftp = Net::FTP->new($host, Debug => 0, Passive => 1) or die "Can't open $host\n";
		#$ftp->login() or die "Cannot login ", $ftp->message;$ftp->cwd($path);$ftp->binary();$ftp->get($blfil,$exe) or die "Failed Blast download: ", $ftp->message;$ftp->quit;
		#sleep(5);
		run_cmd("tar", "-xzf", $exe, "-C", $bdir);
		unlink($exe);
		$exe = "$bdir/ncbi-blast-2.2.29+/bin/blastn";
		@txt = addInfoLtS("blastn",$exe,\@txt,1);
		$exe = "$bdir/ncbi-blast-2.2.29+/bin/makeblastdb";
		@txt = addInfoLtS("makeBlastDB",$exe,\@txt,1);
		}
	}
	if ($installBlast == 2 || $installBlast == 3){
		print "Downloading lambda executables... \n";
		#my $lmdD = "https://lotus2.earlham.ac.uk/lotus/packs/lambda/lambda-v0.9.1-linux_x86-64.tar.gz";
		#if ($isMac){
		#	$lmdD = "https://lotus2.earlham.ac.uk/lotus/packs/lambda/lambda-v0.9.1-darwin_x86-64.tar.gz";
		#}
		if (!$isMac){
			my $lmdD = "https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Linux-x86_64.tar.xz";
			$exe = "$bdir/lambda.tar.xz";
			getS2($lmdD,$exe);
			run_cmd("tar", "-xf", $exe, "-C", $bdir); move("$bdir/lambda3-3.1.0-Linux-x86_64/bin/lambda3", "$bdir/lambda3") or die "Cannot move lambda3: $!\n"; for my $old (bsd_glob("$bdir/lambda3-3*")){ remove_tree($old) if -d $old; unlink($old) if -f $old; }

		}else{
			my $lmdD = "https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Darwin-x86_64.zip";
			$exe = "$bdir/lambda.zip";
			getS2($lmdD,$exe);
			run_cmd("unzip", "-q", "-o", "-d", $bdir, $exe); move("$bdir/lambda3-3.1.0-Darwin-x86_64/bin/lambda3", "$bdir/lambda3") or die "Cannot move lambda3: $!\n"; for my $old (bsd_glob("$bdir/lambda3-3*")){ remove_tree($old) if -d $old; unlink($old) if -f $old; }
		}

		unlink($exe);
		#$exe = "$bdir/lambda/lambda_indexer";
		#@txt = addInfoLtS("lambda_index",$exe,\@txt,1);
		$exe = "$bdir/lambda3";
		@txt = addInfoLtS("lambda3",$exe,\@txt,1);
	}
	#die "$bdir/lambda3";
	if ($installBlast == 0){
		print "\nNo similarity comparison program will be installed.\n";
	}
	#-------BLAST LAMBDA INSTALL END
	
	#usearch
	getUsearch();
	
	#swarm
	print "Downloading swarm executables..\n";
	my $swarmdir = $bdir."swarm-master/";
	my $sexe = "$swarmdir/bin/swarm";
	my $tars = "$bdir/swarm.zip";
	#
	my $swarmtar = "https://lotus2.earlham.ac.uk/lotus/packs/swarm2.1.13.zip";#"https://github.com/torognes/swarm/archive/master.zip";#"https://lotus2.earlham.ac.uk/lotus/packs/swarm206d.tgz";
	getS2($swarmtar,$tars);
	run_cmd("unzip", "-q", "-o", "-d", $bdir, $tars);
	unlink($tars);
	my $callrets = system("make", "-C", "$swarmdir/src/");
	#die($sexe."\n");

	if ($callrets != 0){
		print "\n\n=================\nProblem while compiling swarm.\n"; $finalWarning.="swarm did not compile. The -CL 2 option will not be available to LotuS unless you reinstall swarm manually (lotus.cfg).\n";
	}
	if (-e $sexe){ #not essential
		run_cmd("chmod", "+x", $sexe);
		@txt = addInfoLtS("swarm",$sexe,\@txt,1);
	} else {
		print "Swarm exe did not exist at $sexe\n Therefore swarm was not installed.\n";
	}
	#vsearch: the static Linux x86-64 build shipped as bin/vsearch, otherwise a pinned release
	install_vsearch();

	#infernal
	print "Downloading infernal executables..\n";
	my $iexe = "$bdir/inf112.tar.gz";
	if ($isMac){
		getS2("https://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-macosx-intel.tar.gz",$iexe);
	} else {
		getS2("https://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz",$iexe);
	}
		run_cmd("tar", "-xzf", $iexe, "-C", $bdir);
		$iexe = $isMac ? "$bdir/infernal-1.1.2-macosx-intel/binaries/" : "$bdir/infernal-1.1.2-linux-intel-gcc/binaries/";
	if (-d $iexe){ #not essential
		@txt = addInfoLtS("infernal",$iexe,\@txt,2);
	} else {
		print "infernal binary dir did not exist at $iexe\n Therefore infernal was not installed (fallback to de novo clustal omega).\n";
		$finalWarning .= "infernal binary dir did not exist at $iexe\n Therefore infernal was not installed (fallback to de novo clustal omega).\n";
	}
	unlink("$bdir/inf112.tar.gz") or warn "Could not remove $bdir/inf112.tar.gz: $!\n" if -e "$bdir/inf112.tar.gz";

	#die "$vexe\n";


	#V-Xtractor
	
	my $vxexe = "$bdir/vxtr/vxtractor.pl";
	ensure_dir("$bdir/vxtr/");
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/VXtractor/vxtractor.pl",$vxexe);
	@txt = addInfoLtS("vxtractor",$vxexe,\@txt,1);
	$vxexe = "$bdir/vxtr/HMM.zip";
	getS2("https://lotus2.earlham.ac.uk/lotus/packs/VXtractor/HMMs.zip",$vxexe);
	print("unzip -o -q $vxexe -d $bdir/vxtr/;rm $vxexe;");
	run_cmd("unzip", "-o", "-q", $vxexe, "-d", "$bdir/vxtr/"); unlink($vxexe) or warn "Could not remove $vxexe: $!\n";
	@txt = addInfoLtS("vxtractorHMMs","$bdir/vxtr/HMMs/",\@txt,2);
	#die "$bdir/vxtr/HMMs/";




	## iqtree2
	print "Downloading IQ-TREE 2 executables..\n";
	$dtar = "$bdir/iqtree-2.1.1-Linux.tar.gz";
	$dexe = "$bdir/iqtree-2.1.1-Linux/bin/iqtree2";
	if ($isMac){
		getS2("https://github.com/iqtree/iqtree2/releases/download/v2.1.1/iqtree-2.1.1-MacOSX.zip",$dtar);
		$dexe = "$bdir/iqtree-2.1.1-MacOSX/bin/iqtree2";
	} else {
		getS2("https://github.com/iqtree/iqtree2/releases/download/v2.1.1/iqtree-2.1.1-Linux.tar.gz",$dtar);
	}
	if ($isMac){ run_cmd("unzip", "-q", "-o", "-d", $bdir, $dtar); } else { run_cmd("tar", "-xzf", $dtar, "-C", $bdir); }
	unlink($dtar);

	if (-e $dexe){ #not essential
		run_cmd("chmod", "+x", $dexe);
		@txt = addInfoLtS("iqtree",$dexe,\@txt,1);
	} else {
		$finalWarning .= "iqtree2 exe did not exist at $dexe\n Therefore iqtree2 was not installed (please manually install).\n";
		print "iqtree2 exe did not exist at $dexe\n Therefore iqtree2 was not installed (please manually install).\n";
	}

	##mafft
	print "Downloading MAFFT 7 executables..\n";
	$dtar = "$bdir/mafft-7.471-linux.tgz";
	$dexe = "$bdir/mafft-linux64/mafft.bat";
	if ($isMac){
		getS2("https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip",$dtar);
		$dexe = "$bdir/mafft-mac/mafft.bat";
		run_cmd("unzip", "-q", "-o", "-d", $bdir, $dtar);

	} else {
		getS2("https://mafft.cbrc.jp/alignment/software/mafft-7.471-linux.tgz",$dtar);
		run_cmd("tar", "-xzf", $dtar, "-C", $bdir);
	}
	unlink($dtar);

	if (-e $dexe){ #not essential
		run_cmd("chmod", "+x", $dexe);
		@txt = addInfoLtS("mafft",$dexe,\@txt,1);
	} else {
		$finalWarning .= "MAFFT exe did not exist at $dexe\n Therefore MAFFT was not installed (please manually install).\n";
		print "MAFFT exe did not exist at $dexe\n Therefore MAFFT was not installed (please manually install).\n";
	}

	#fasttree
	print "Downloading FastTree executables..\n";
	$exe = "$bdir/FastTreeMP";
	my $exe1 = "$bdir/FastTree.c";
	#system("wget -O $exe $fastt");
	#my $fastt = "http://www.microbesonline.org/fasttree/FastTreeMP";
	#if ($isMac){}
	my $fastt = "https://lotus2.earlham.ac.uk/lotus/packs/FastTree.c"; #http://www.microbesonline.org/fasttree/
	getS2($fastt,$exe1);
	my $cc = command_exists("gcc") // command_exists("cc") // "gcc";
	$callret = system($cc, "-DOPENMP", "-fopenmp", "-O3", "-finline-functions", "-funroll-loops", "-Wall", "-o", $exe, $exe1, "-lm");
	if ($callret != 0){
		print "\n\n=================\nProblem while compiling fasttree, trying fasttree without multithread and SSE support (might be slower, but if it's working..)\n";
		$finalWarning .= "fasttree compiled without multithreading support (you can not use the -thr LotuS option.\n";
		$exe = "$bdir/FastTree";
		$callret = system($cc, "-DNO_SSE", "-O3", "-finline-functions", "-funroll-loops", "-Wall", "-o", $exe, $exe1, "-lm");}
	if ($callret != 0){
		#the remaining programs still install; only FastTree trees (-buildPhylo 1, the default) are unavailable
		my $msg = "fasttree compilation failed, so FastTree was not installed. This is most likely an issue with your C compiler or the OpenMP libraries (see http://www.microbesonline.org/fasttree/#Install). Until it is installed, run LotuS3 with -buildPhylo 2 (IQ-TREE) or -buildPhylo 0.\n";
		print "\n\n=================\n$msg"; $finalWarning .= $msg;
	} else {
		run_cmd("chmod", "+x", $exe);
		@txt = addInfoLtS("fasttree",$exe,\@txt,1);
	}


	#flash
	if (0){#not needed in lotus2 any longer..
		my $flashdir = $bdir."FLASH-1.2.10";
		my $fexe = "$flashdir/flash";
		my $tar = "$bdir/Flash.tar.gz";
		my $flashTar = "https://lotus2.earlham.ac.uk/lotus/packs/FLASH-1.2.10.tar.gz";#"http://sourceforge.net/projects/flashpage/files/FLASH-1.2.10.tar.gz/download";
		getS2($flashTar,$tar);
		run_cmd("tar", "-xzf", $tar, "-C", $bdir);
		unlink($tar);
		$callret = system("make", "-C", $flashdir);
		if ($callret != 0){
			print "\n\n=================\nProblem while compiling FLASH.\n"; $finalWarning.="Flash did not compile. This means you can not use paired reads with LotuS.\n";
		}
		run_cmd("chmod", "+x", $fexe);
		@txt = addInfoLtS("flashBin",$fexe,\@txt,1);
	}
	
	


	#cd-hit
	my $cdhitdir = $bdir."cdhit-master/";
	my $cexe = "$cdhitdir/cd-hit-est";
	my $ctar = "$bdir/cdhit.zip";
	#my $cdhitTar = "https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz";
	my $cdhitTar = "https://lotus2.earlham.ac.uk/lotus/packs/cd-hit_git.zip";#"https://github.com/weizhongli/cdhit/archive/master.zip";
	getS2($cdhitTar,$ctar);
	#system("tar -xzf $tar -C $bdir");
	run_cmd("unzip", "-o", "-q", $ctar, "-d", $bdir);
	unlink($ctar);
	$callret = system("make", "-C", $cdhitdir);
	if ($callret != 0){
		print "\n\n=================\nProblem while compiling CD-HIT.\n"; $finalWarning.="CD-HIT did not compile. The -UP 3 option will not be available to LotuS unless you reinstall cd-hit-est manually (and add to lotus.cfg). \n";
	} else {
		run_cmd("chmod", "+x", $cexe);
		@txt = addInfoLtS("cd-hit",$cexe,\@txt,1);
	}


	my $rdpf = "https://lotus2.earlham.ac.uk/lotus/packs/rdp_classifier_2.12.zip"; #"http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.6.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Frdp-classifier%2F&ts=1391590725&use_mirror=netcologne";
	#RDP classifier
	$exe = "$bdir/rdp.zip";
	#system("wget -O $exe $rdpf");
	getS2($rdpf,$exe);
	#die("unzip $exe -d $bdir");
	run_cmd("unzip", "-o", "-q", $exe, "-d", $bdir);
	unlink($exe);
	$exe = $bdir."rdp_classifier_2.12/dist/classifier.jar";
	@txt = addInfoLtS("RDPjar",$exe,\@txt,1);



	#clustalO
	#optional: LotuS aligns with MAFFT. The former macOS binary is gone upstream (HTTP 403).
	if ($isMac){
		print "Clustal Omega was not installed on macOS (no verified binary available); it is optional, MAFFT is used for alignments.\n";
		return;
	}
	my $clo = "https://lotus2.earlham.ac.uk/lotus/packs/clustalo-1.2.0-Ubuntu-x86_64";
	$exe = "$bdir/clustalo-1.2.0-Ubuntu-x86_64";
	getS2($clo,$exe);
	run_cmd("chmod", "+x", $exe);
	@txt = addInfoLtS("clustalo",$exe,\@txt,1);
}

sub user_options(){

	if ($condaDBinstall){#no user input at all wanted
		return;
	}
	if ( $UID ne "??" || $usearchInstall ne ""){#a configured UID is sufficient to identify a previous installation
		my $inp="";
		
		if ($usearchInstall eq ""){
			while ($inp !~ m/^[123]$/){
				print "Detected previous installation of LotuS, do you want to \n";
				#print " (1) search & install updates\n";
				print " (1) refresh databases and reinstall secondary software (e.g. after \"git pull\")\n";
				print " (2) refresh only databases (secondary software remains unchanged)\n";
				print " (3) set or update the path to your USEARCH binary\n";
				print "Answer: \n";
				$inp = read_user_input("the previous-installation choice");
			}
		}
		if ($inp eq "3"){
			print "Enter the full (absolute) path to your usearch binary:\n";
			while ($usearchInstall eq ""){
				$usearchInstall = read_user_input("the USEARCH path");
				if (!-f $usearchInstall || !-x $usearchInstall){
					$usearchInstall="";
					print "The path is not an executable file; please re-enter it (or abort with Ctrl-c):\n";
				} else {
					$usearchInstall = abs_path($usearchInstall);
				}
			}
		}
		if ($usearchInstall ne ""){
			print "Setting usearch binary (required for lotus) to \n";
			if (!-f $usearchInstall || !-x $usearchInstall){die "USEARCH path $usearchInstall is not an executable file.\n";}
			$usearchInstall = abs_path($usearchInstall);
			@txt = addInfoLtS("usearch",$usearchInstall,\@txt,1);
			print "Successfully added usearch into LotuS. Now LotuS is ready to run.\n";
			finishAI("none");
			exit(0);
		}
		if ($inp eq "2"){
			$onlyDbinstall = 1;
		}
	}
	#auto update END
	my $installAll = 0;
	if ($onlyDbinstall){
		print "Installing LotuS tax databases anew.. \nplease choose which databases to install in the following dialogs\n\n";
	} else {
		print "Some programs require a recent C++ compiler. Existing files are retained until their replacements download successfully, and lOTUs.cfg will be updated.\n";
		print "Install LotuS3 with all possible dependencies (all databases, ITS, ONT related workflows)?\nSimply enter or \"1\" for yes, \"0\" for detailed configuration via question.\nAnswer: ";
		$installAll = read_user_input("the all-dependencies choice", [0, 1], 1);
		if ($isMac){print "Mac system detected, installing corresponding mac software.\n";}
		if ($installAll) {
			$installBlast = 3; #both BLAST and Lambda
			@refDBinstall = (0) x 10;
			$refDBinstall[8] = 1; #all supported similarity-reference databases
			$ITSready = 1;
			$getUTAX = 1;
			$installONT = 1;
			print "Selected all supported dependencies: all reference databases, ITS/UTAX resources, BLAST and Lambda, ONT tools, and the standard programs/R packages.\n";
		} else {
			print "\n\nFor similarity based taxonomic assignments LotuS can either use \n (1) Blastn \n (2) Lambda \n (3) both, decide at runtime which to use or\n (0) none\n Answer:";
			$installBlast = read_user_input("the similarity-search program choice", [0, 1, 2, 3]);
		}
	}

	if (!$installAll) {
		print "\n\nDo you want to install a reference database for similarity based annotations?\n";
		print " (1) KSGP (~1.5 GB), covering SSU for Archaea, Bacteria and Eukaryotes, 2026 release. \n (2) SILVA (~2.5 GB), contains LSU as well as SSU, 138.1 2020 release.\n (3) GreenGenes2 (~1 GB), 2022 release.\n (4) HITdb (~100 MB) 16S bacterial database specialized on the gut environment.\n";
		print " (5) PR2 (~100 MB), an SSU database specialized for marine eukaryotes.\n";
		print " (6) beeTax (~2 MB) database specialized (and named) on taxonomy specific to the bee gut.\n";
		print " (8) KSGP + SILVA + GG2 + PR2 + HITdb + beeTax (select a specific DB in each LotuS3 run)\n (0) no database.\n";
		print "Answer:";
		my $choice = read_user_input("the reference-database choice", [0 .. 6, 8]);
		@refDBinstall = (0) x 10;
		$refDBinstall[$choice] = 1;
	}

	#SILVA license
	if ($refDBinstall[2] || $refDBinstall[8]){
		print "Please read the SILVA license: https://www.arb-silva.de/fileadmin/silva_databases/LICENSE.txt. Do you accept (y/n)? \n";
		while (1){
			my $choice = lc(read_user_input("the SILVA license response"));
			if ($choice eq "y" || $choice eq "yes"){
				last;
			} elsif ($choice eq "n" || $choice eq "no") {
				die "You need to accept the SILVA license before installation can continue.\n";
			}
			print "Please answer y or n: ";
		}
	}

	return if $installAll; #the license response above is still required

	print "\n\n -- ITS -- Do you want to\n (1) install databases and programs required to process ITS data (including fungi ITS UNITE database)\n (0) no ITS related packages\n Answer:";
	$ITSready = read_user_input("the ITS package choice", [0, 1]);

	print "\n\n -- UTAX -- Do you want to\n (1) install utax taxonomic classification databases (16S, ITS)?\n (0) no utax related databases\n Answer:";
	$getUTAX = read_user_input("the UTAX database choice", [0, 1]);

	if (!$onlyDbinstall) {
		print "\n\n -- ONT -- Install ONT tools (Savont and Barbell)?\nSimply enter or \"1\" for yes, \"0\" for no.\nMinimap2 is included in the core installation for read mapping.\nAnswer: ";
		$installONT = read_user_input("the ONT tools choice", [0, 1], 1);
	}
}
