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
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
sub addInfoLtS;sub finishAI;
#subroutines to download various DBs..
sub getGG; sub getGG2; sub getSLV;sub getHITdb; sub getPR2db;sub getKSGP;sub getbeetax;
sub buildIndex;
sub get_DBs;
sub getS2;
sub getInfoLtS;
sub getInstallVer;
sub compile_sdm;
sub compile_LCA;
sub compile_rtk;
sub checkLtsVer;
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
my $FILEfetch = eval{
  require File::Fetch;
  File::Fetch->import();
  1;
};
my $LWPsimple = eval{
  require LWP::Simple;
  LWP::Simple->import();
  1;
};
my $forceUpdate=0;
my $condaDBinstall=0;
my $downloadLmbdIdx = 0; #download lambda index from webpage
my $compile_lambda=0;
my $usearchInstall = "";
my $noTelemetry = 0;

GetOptions(
	"forceUpdate"     => \$forceUpdate,
	"condaDBinstall"  => \$condaDBinstall,
	"downloadLmbdIdx" => \$downloadLmbdIdx,
	"lambdaIndex"     => \$compile_lambda,
	"link_usearch=s"  => \$usearchInstall,
	"no-telemetry"    => \$noTelemetry,
) or die "Invalid command line options\n";

if ($compile_lambda && $downloadLmbdIdx){
	die "Can't use both -lambdaIndex and -downloadLmbdIdx arguments together\nAborting..\n";
}

my $WGETpres = command_exists("wget") ? 1 : 0;

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
my $UID = getInfoLtS("UID",\@txt);
my $uspath = getInfoLtS("usearch",\@txt);

#usearch binary linking is handled by GetOptions above


##### TESTING / DEBUG ##########
#		@txt = getKSGP(\@txt);die;#@txt = getGG2(\@txt); 

#DEBUG
#@txt = getPR2db(\@txt);;exit;

###### GET USER OPTIONS ###################

my ($lver,$sver) = getInstallVer("");
if ($forceUpdate==0 && $condaDBinstall == 0){
	print "\n\t####################################\n\t LotuS $lver Auto Installer script.\n\t####################################\n\n";
} elsif ($condaDBinstall){
	print "\n\nConda LotuS install: downloading the standard database set for LotuS3\n\n";
} else {
	print "\n\nRerunning updates due to updated autoupdate.pl script\n\n";
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


###################   database downloads ... #########################
get_DBs();

if ($condaDBinstall){
	finishAI("d");
	print "Finished LotuS3 database install (Conda autoinstall)\nEnjoy LotuS3!\n";
	exit(0);
}

if ($onlyDbinstall){
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
#first install/compile sdm,LCA,rtk
my $nsdmp = compile_sdm("$ldir/sdm_src");
@txt = addInfoLtS("sdm",$nsdmp,\@txt,1);
$nsdmp = compile_LCA("$ldir/LCA_src");
@txt = addInfoLtS("LCA",$nsdmp,\@txt,1);

$nsdmp = compile_rtk("$ldir/rtk_src");
@txt = addInfoLtS("rtk",$nsdmp,\@txt,1);

#download and install the remaining programs exactly once
get_programs();

finishAI("");


print "\n\nInstallation script finished.\nPlease read the README for examples and references to proprietary software used in this pipeline.\n";


#After install on your system, open\n   ".$ldir."lOTUs.cfg\nand search for the entry \"usearch {XX}\".\nReplace {XX} with the absolute path to your usearch install, e.g. /User/Thomas/bin/usearch/usearch7.0.1001_i86linux32\n LotuS is ready to run.\n";

sub finishAI($){
	my ($vTag) = @_;
	#write new cfg file
	write_config_atomic("$ldir/lOTUs.cfg", \@txt);
	return if ($vTag eq "none");
	if ($LWPsimple && !$noTelemetry){
		my $external_php = get("https://lotus2.earlham.ac.uk/lotus/in.php?ID=$UID&VERSION=$vTag$lver&SDMV=$sver") || print "";
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
		my $sdmF = "$sdmsrc/IO.h";
		if (-e $sdmF){
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
	while ($txt[$i] !~ m/^$ss\s/){
		#print $txt[$i]."\n";
		$i++;
		if ($i >= @txt){
			#die ("Could not find the entry \"$cmd\" in lotus configuration file. Aborting Installer..\n")
			print "Could not find the entry \"$cmd\" in lotus configuration file. Inserting anew..\n";
			push(@txt,""); last;
		}
	}
	$txt[$i] = $cmd." ".$ex."\n";
	$i++;
	while ($i<@txt){ if ($txt[$i] =~ m/^$ss\s/){splice(@txt,$i,1) ; $i--;} $i++; last if ($i>=@txt); }
	print "done.\n";
	
	write_config_atomic("$ldir/lOTUs.cfg", \@txt);

	#DEBUG
	#print $txt[$i]."\n";
	return @txt;
}
               
sub getInfoLtS($ $){
	my ($cmd,$aref) = @_;
	my $ss = quotemeta $cmd;
	foreach my $line (@{$aref}){
		chomp(my $copy = $line);
		return $1 if ($copy =~ m/^$ss\s+(.*)$/);
		return "??" if ($copy =~ m/^$ss\s*$/);
	}
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
	open I,"<$DBin" or die "Can;t open PR2 fasta $DBin\n";
	open F,">$DBin.tmp" or die "Can;t open PR2 fasta tmp $DBin.tmp\n";
	#>AB353770.1.1740_U;tax=k:Eukaryota,d:TSAR,p:Alveolata-Dinoflagellata,c:Dinophyceae,o:Peridiniales,f:Kryptoperidiniaceae,g:Unruhdinium,s:Unruhdinium_kevei

	while (my $l = <I>){
		chomp $l;
		if ($l =~ m/^>/){
			 my @spl = split /;tax=/,$l;
			$spl[0] =~ s/^>//;
			print F ">".$spl[0]."\n";
			#my @spl2 = split(/;/,$spl[1]);
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
	rename($taxTmp, $tout) or die "Cannot replace $tout with $taxTmp: $!\n";
	unlink($DBin) if -e $DBin; move("$DBin.tmp", $DBin) or die "Cannot replace $DBin with $DBin.tmp: $!\n";
}


sub buildIndex($){
	my ($DBfna) = @_;
	return unless ($compile_lambda);
	
	my $lambdaIdxBin = "";#find where lambda is installed in
	foreach my $line (@txt){
		if ($line =~ m/^lambda3\s+(\S+)/ ) {$lambdaIdxBin = $1;}
	}
	if ($lambdaIdxBin ne "" && !-x $lambdaIdxBin){
		$lambdaIdxBin = command_exists($lambdaIdxBin) // "";
	}
	die "Cannot build Lambda index: database file $DBfna is missing or empty\n" unless (-s $DBfna);
	die "Cannot build Lambda index: lambda3 is not configured or executable\n" unless ($lambdaIdxBin ne "" && -x $lambdaIdxBin);
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
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.fna",$DB);
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/beeTax_Engel/beEngel.txt",$DBtax);
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
	#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/gb203PR2.tar.gz",$DB.".tar.gz");
#	system "tar -xzf $DB.tar.gz -C $ddir/PR2;rm $DB.tar.gz";
	#getS2("https://lotus2.earlham.ac.uk/lambdaDBs/v3.0/gb203_pr2_all_10_28_99p.fasta.lba.gz","$ddir/PR2/gb203_pr2_all_10_28_99p.fasta.lba.gz") if ($downloadLmbdIdx);

#https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_mothur.tax.gz

	#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/PR2/pr2_version_5.0.0_SSU_mothur.tax.gz",$DBtax.".gz");
	
	#my $DB = "$ddir/PR2_5.0_pre.fasta";
	my $DB = "$ddir/PR2_5.0.fasta";
	#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/PR2//pr2_version_5.0.0_SSU_mothur.fasta.gz",$DB.".gz");
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/PR2//pr2_version_5.0.0_SSU_UTAX.fasta.gz",$DB.".gz");
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
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_sequences.fna",$DB);
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/hitdb/HITdb_taxonomy_qiime.txt",$DBtax.".pre");
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


sub getGG($){
	my ($aref) = @_;
	die "getGG::Greengenes is no longer supported\n";
	my @txt = @{$aref};
	#greengenes ------------------------
	my $gg1 = "http://lotus2.earlham.ac.uk/lotus/packs/gg_13_5.fasta.gz";
	my $gg2 = "http://lotus2.earlham.ac.uk/lotus/packs/gg_13_5_taxonomy.gz";
	my $DB = "$ddir/gg_13_5.fasta";
	unlink glob("${DB}*");
	#system("wget -O $DB.gz $gg1");
	print "Downloading Greengenes may 2013 release..\n";
	getS2($gg1,"$DB.gz");
	sleep(10);
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("TAX_REFDB_GG",$DB,\@txt,1);
	buildIndex($DB);
	$DB = "$ddir/gg_13_5_taxonomy";
	#system("wget -O $DB.gz $gg2");
	getS2($gg2,"$DB.gz");
	sleep(3);
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("TAX_RANK_GG",$DB,\@txt,1);
	return @txt;
}

sub getKSGP($){
	my ($aref) = @_;
	my @txt = @{$aref};
	

	my $DB = "$ddir/KSGPv4.0";
	print "Downloading KSGP v4.0 Jul 2026 release..\n";
	my $tarUTN = "$ddir/KSGPv4.0.gz";	my $tarUTNtax = "$ddir/KSGPv4.0.tax.gz";
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_v4.0.fasta.gz",$tarUTN);
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGP_plus2.tax.gz",$tarUTNtax);
	gunzip_file($tarUTN, "$DB.fasta"); gunzip_file($tarUTNtax, "$DB.tax");
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv4.0/KSGPv4.0.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);



#	my $DB = "$ddir/KSGP_v3.1";unlink glob("${DB}*");
#	print "Downloading KSGP v3.1 2025 release..\n";
#	my $tarUTN = "$ddir/KSGPv3.1.gz";	my $tarUTNtax = "$ddir/KSGPv3.1.tax.gz";
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1/KSGP.fasta.gz",$tarUTN);
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1/KSGP.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3.1//KSGP_v3.1.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);

#	my $DB = "$ddir/KSGP_v2.0";unlink glob("${DB}*"); print "Downloading KSGP v3 2025 release..\n";
#	my $tarUTN = "$ddir/KSGPv3.gz";	my $tarUTNtax = "$ddir/KSGPv3.tax.gz";
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3/KSGP_v3.fasta.gz",$tarUTN);
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv3/KSGP_v3.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");


#	print "Downloading KSGP 2024 release..\n";
#	my $DB = "$ddir/KSGP_v2.0";unlink glob("${DB}*"); my $tarUTN = "$ddir/KSGPv2.gz";	my $tarUTNtax = "$ddir/KSGPv2.tax.gz";
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2/KSGP_v2.fasta.gz",$tarUTN);
#	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2/KSGP_LCA_v2.tax.gz",$tarUTNtax);
#	system("gunzip -c $tarUTN > $DB.fasta");system("gunzip -c $tarUTNtax > $DB.tax");
#	system("rm -f $tarUTN $tarUTNtax");
	
	
	#getS2("https://ksgp.earlham.ac.uk/downloads/v1.0/KSGP_v1.0.tar.gz",$tarUTN);
	#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/KSGPv2//KSGP_v2.0.fasta.lba.gz","$DB.fasta.lba.gz") if ($downloadLmbdIdx);
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
		$baseSP = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/";
		#my $SlvAltFna = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_SSU.fasta.gz";
		#getS2($SlvAltFna,"$DB.gz");
		#system("gunzip -c $DB.gz > $DB;rm -f $DB.gz"); 
		#my $SlvAltTax = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_SSU.tax.gz";
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
		$baseSP = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/";
	}
	#	my $SlvAltFna = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_LSU.fasta.gz";
	#	getS2($SlvAltFna,"$DB.gz");
	#	system("gunzip -c $DB.gz > $DB;rm -f $DB.gz"); 
	#	my $SlvAltTax = "http://lotus2.earlham.ac.uk/lotus/packs/DB/SLV/SLV_132_LSU.tax.gz";
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


sub getS2($ $){
	my ($in,$out) = @_;
	print "getS2:$in\n$out\n";
	die "Refusing non-http(s) download URL: $in\n" unless ($in =~ m{^https?://}i);
	ensure_dir(dirname($out));
	my $tmp = "$out.tmp.$$";
	unlink($tmp) if (-e $tmp);
	if ($WGETpres){
		print "wget -O $tmp $in\n";
		run_cmd("wget", "-O", $tmp, $in);
	} elsif (!$isMac && $LWPsimple){
		print "LWP\n";
		my $rc = getstore($in,$tmp);
		die "Download failed for $in: HTTP/status $rc\n" unless (defined($rc) && $rc >= 200 && $rc < 300);
	} elsif ($FILEfetch){
		print "FETCH\n";
		my $ff = File::Fetch->new( uri => $in);
		my $file = $ff->fetch() or die "Can't download file $in with File::Fetch\n".$ff->error()."\n";
		move($file, $tmp) or die "Can't move fetched file $file to $tmp: $!\n";
	} else {
		die "no suitable library / program on you system. Please ensure that \"wget\" is installed\n";
	}
	die "Download produced no file: $in -> $tmp\n" unless (-e $tmp);
	die "Downloaded file is empty: $in -> $tmp\n" unless (-s $tmp);
	rename($tmp, $out) or do {
		unlink($tmp) if (-e $tmp);
		die "Can't replace $out with downloaded file $tmp: $!\n";
	};
	return $out;
}

sub checkLtsVer($){
	my ($lver) = @_;
	die "LWP::Simple is required for the updater\n" if (!$LWPsimple);
	my $updtmpf = get("http://lotus2.earlham.ac.uk/lotus/lotus/updates/Msg.txt");
	die "Could not download update message list\n" unless defined($updtmpf);
	my $msg = ""; my $hadMsg=0;
	open( TF, '<', \$updtmpf ); while(<TF>){$msg .= $_;}  close(TF); 
	foreach my $lin (split(/\n/,$msg)){
		my @spl = split /\t/,$lin;
		next if (@spl==0);
		if (version_is_newer($spl[0],$lver)){print $spl[1]."\n\n"};
		$hadMsg=1;
	}
	# compare to server version
	$updtmpf = get("http://lotus2.earlham.ac.uk/lotus/lotus/updates/curVer.txt");
	die "Could not download current LotuS version\n" unless defined($updtmpf);
	open( TF, '<', \$updtmpf ); my $lsv = <TF>; close(TF); chomp $lsv;
	die "Updater returned an invalid version '$lsv'\n" unless ($lsv =~ m/^\d+(?:\.\d+)+$/);
	my $msgEnd = "";
	$updtmpf = get("http://lotus2.earlham.ac.uk/lotus/lotus/updates/curVerMsg.txt");
	die "Could not download current update message\n" unless defined($updtmpf);
	open( TF, '<', \$updtmpf ); while(<TF>){$msgEnd .= $_;} close(TF); 
	
	$updtmpf = get("http://lotus2.earlham.ac.uk/lotus/lotus/updates/UpdateHist.txt");
	die "Could not download update history\n" unless defined($updtmpf);
	my $updates = "";
	open( TF, '<', \$updtmpf );$msg = ""; while(<TF>){$msg .= $_;}  close(TF); 
	foreach my $lin (split(/\n/,$msg)){
		my @spl = split /\t/,$lin; chomp $lin;
		next if (@spl < 2 || $spl[0] eq "");
		if ($spl[1] =~ m/LotuS (\d+(?:\.\d+)+)/){
			if (version_is_newer($1,$lver)){$updates.= $spl[0]."\t".$spl[1]."\n"};
		}
	}
	if ($updates ne ""){
		print "--------------------------------\nThe following updates are available:\n--------------------------------\n";
		print $updates;
		print "\n\nCurrent Lotus version is :$lver\nLatest version is: $lsv\n";
	}
	
	if ($hadMsg || $updates ne ""){sleep(4);}

	
	#die;
	return $lsv,$msgEnd;
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
		unlink glob("$ldi2/*.o");
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
sub compile_rtk($){
	my ($ldi2) = @_;
	my $expPath = "$bdir/rtk";
	if (-x $expPath){#test if can execute locally
		my ($rtkV,$status) = capture_cmd($expPath, "-v");
		return $expPath if ($status == 0 && $rtkV =~ m/rtk \d/);
	}

	if (-d $ldi2 && -f "$ldi2/Makefile" ){
		print "Compiling rtk..\n";
		unlink glob("$ldi2/*.o");
		my $stat = system("make", "-C", $ldi2);
		if ($stat == 0){
			unlink("$ldir/rtk") if -e "$ldir/rtk"; unlink("$bdir/rtk") if -e "$bdir/rtk"; move("$ldi2/rtk", "$bdir/rtk") or die "Cannot install rtk: $!\n"; run_cmd("chmod", "+x", "$bdir/rtk");
		} else {
			die "Compilation of required rtk binary failed (make status $stat). Install a C++ compiler and rerun the installer.\n";
		}
	} else {
		die "rtk source directory or Makefile is missing at $ldi2\n";
	}
	die "Compilation did not produce executable $expPath\n" unless (-e $expPath);
	run_cmd("chmod", "+x", $expPath);
	return $expPath;
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
		unlink glob("$ldi2/*.o");
		my $stat = system("make", "-C", $ldi2);
		if ($stat != 0){#repeat without gzip
			print "\n\n\n\n=================\nProblem compiling sdm with gzip support\nFallback to sdm compilation without gzip support\n";
			my $header = "$ldi2/DNAconsts.h";
			my $backup = "$header.installer-backup.$$";
			copy($header, $backup) or die "Cannot back up $header before fallback compilation: $!\n";
			run_cmd($^X, "-pi", "-e", "s/#define _gzipread/#define _notgzip/g", $header);
			unlink glob("$ldi2/*.o");
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
	my $status = $raw_status == -1 ? -1 : ($raw_status >> 8);
	return ($output,$status);
}

sub read_user_input {
	my ($context) = @_;
	my $line = <STDIN>;
	die "End of input while waiting for $context; installation aborted.\n" unless defined($line);
	chomp($line);
	return $line;
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
	my ($check,$status) = capture_cmd($exe, "--version");
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
		getS2("http://drive5.com/utax/data/utax_rdp_16s_tainset15.tar.gz",$tarUTN);
		run_cmd("tar", "-xzf", $tarUTN, "-C", $ddir); unlink($tarUTN) or warn "Could not remove $tarUTN: $!\n";
		$tarUTN="$ddir/utax_ITS.tar.gz";
		getS2("http://drive5.com/utax/data/utax_unite_v7.tar.gz",$tarUTN);
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
		#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_refs_qiime_ver8_99_s_all_02.02.2019.fasta.zip",$tarUN);
#		my $UNITEdb = "$ddir/UNITE/sh_refs_qiime_ver8_99_s_all_02.02.2019.fasta";
		#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/sh_qiime_release_02.03.2015.zip",$tarUN);
		#system("rm -fr $ddir/UNITE;unzip -q -o $tarUN -d $ddir/UNITE/");
		#getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/sh_taxonomy_qiime_ver8_99_s_all_02.02.2019.txt.zip",$tarUN);
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
		getS2("https://lotus2.earlham.ac.uk/packs/DB/UNITE/Lambda3/sh_refs_v10_19.02.2025.fasta.lba.gz","$UNITEdb.fasta.lba.gz") if ($downloadLmbdIdx);
		buildIndex("$UNITEdb.fasta");
		
		@txt = addInfoLtS("TAX_REFDB_ITS_UNITE","$UNITEdb.fasta",\@txt,1);
		@txt = addInfoLtS("TAX_RANK_ITS_UNITE","$UNITEdb.tax",\@txt,1);

		unlink($tarUN);
	}
	
	if ($ITSready){#ITS chimera check ref DB
		#my $itsDB = "http://lotus2.earlham.ac.uk/lotus/packs/DB/uchime_reference_dataset_11.03.2015.zip";
		my $itsDB = "http://lotus2.earlham.ac.uk/lotus/packs/DB/UNITE/uchime/uchime_UNITE_16_10_22.zip";
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
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/DB/phiX.fasta",$phiXf);
	@txt = addInfoLtS("REFDB_PHIX",$phiXf,\@txt,1);
	
	#db gold #exchanged for rdp_gold since 1.30
	#my $goldDB = "http://drive5.com/uchime/gold.fa";
	my $goldDB = "http://lotus2.earlham.ac.uk/lotus/packs/rdp_gold.fa.gz";
	my $DB = "$ddir/rdp_gold.fa";
	#system("wget -O $DB $goldDB");
	getS2($goldDB,$DB.".gz");
	gunzip_file("$DB.gz", $DB);
	@txt = addInfoLtS("UCHIME_REFDB",$DB,\@txt,1);



	#db Silva 119 clustered to 93% for LSUs
	my $LTUrefDB = "http://lotus2.earlham.ac.uk/lotus/packs/SILVA_119_LSU_93.ref.fasta.gz";
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


sub get_programs{
	#-----------  exit prog here, if set
	#-----------------------



	#minimap2
	print "Installing minimap2 executable..\n";
	my $dtar = "$bdir/minimap2-2.28_x64-linux.tar.bz2";
	my $dexe = "$bdir/minimap2-2.28_x64-linux/minimap2";
	if ($isMac){
		$dexe = command_exists("minimap2") // "";
		if ($dexe ne ""){
			@txt = addInfoLtS("minimap2",$dexe,\@txt,1);
		} else {
			my $message = "minimap2 was not found in PATH. Install it natively on macOS and rerun the installer.\n";
			print $message;
			$finalWarning .= $message;
		}
	} else {
		getS2("http://lotus2.earlham.ac.uk/lotus/packs//minimap2-2.28_x64-linux.tar.bz2",$dtar);
		run_cmd("tar", "-xjf", $dtar, "-C", $bdir);
		unlink($dtar) or warn "Could not remove $dtar: $!\n";
		if (-e $dexe){ #not essential
			run_cmd("chmod", "+x", $dexe);
			@txt = addInfoLtS("minimap2",$dexe,\@txt,1);
		} else {
			$finalWarning .= "minimap2 executable did not exist at $dexe; minimap2 was not installed.\n";
			print "minimap2 executable did not exist at $dexe; minimap2 was not installed.\n";
		}
	}


	if ($ITSready){ #ITSx
		#itsx
		print "Downloading ITSX to detect valid ITS regions..\n";
		my $tarUTN = "$bdir/ITSx_1.1.4.tar.gz";
		getS2("http://lotus2.earlham.ac.uk/lotus/packs/ITSx_1.1.4.tar.gz",$tarUTN);
		run_cmd("tar", "-xzf", $tarUTN, "-C", $bdir); unlink($tarUTN) or warn "Could not remove $tarUTN: $!\n";
		@txt = addInfoLtS("itsx","$bdir/ITSx_1.1.4/./ITSx",\@txt,1);
		@txt = addInfoLtS("hmmsearch","$bdir/ITSx_1.1.4/bin/hmmsearch",\@txt,1);

	}

	#-------BLAST LAMBDA INSTALL
	if ($installBlast == 1 || $installBlast == 3){
		#Blast
		print "Downloading blast executables...\n";
		my $blfil = "ncbi-blast-2.2.29+-x64-linux.tar.gz";
		if ($isMac){
			$blfil = "ncbi-blast-2.2.29+-universal-macosx.tar.gz";
		}

		$exe = "$bdir/blast.tar.gz";
		getS2("http://lotus2.earlham.ac.uk/lotus/packs/".$blfil,$exe);
			
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
	if ($installBlast == 2 || $installBlast == 3){
		print "Downloading lambda executables... \n";
		#my $lmdD = "http://lotus2.earlham.ac.uk/lotus/packs/lambda/lambda-v0.9.1-linux_x86-64.tar.gz";
		#if ($isMac){
		#	$lmdD = "http://lotus2.earlham.ac.uk/lotus/packs/lambda/lambda-v0.9.1-darwin_x86-64.tar.gz";
		#}
		if (!$isMac){
			my $lmdD = "https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Linux-x86_64.tar.xz";
			$exe = "$bdir/lambda.tar.xz";
			getS2($lmdD,$exe);
			run_cmd("tar", "-xf", $exe, "-C", $bdir); move("$bdir/lambda3-3.1.0-Linux-x86_64/bin/lambda3", "$bdir/lambda3") or die "Cannot move lambda3: $!\n"; for my $old (glob("$bdir/lambda3-3*")){ remove_tree($old) if -d $old; unlink($old) if -f $old; }

		}else{
			my $lmdD = "https://github.com/seqan/lambda/releases/download/lambda-v3.1.0/lambda3-3.1.0-Darwin-x86_64.zip";
			$exe = "$bdir/lambda.zip";
			getS2($lmdD,$exe);
			run_cmd("unzip", "-q", "-o", "-d", $bdir, $exe); move("$bdir/lambda3-3.1.0-Darwin-x86_64/bin/lambda3", "$bdir/lambda3") or die "Cannot move lambda3: $!\n"; for my $old (glob("$bdir/lambda3-3*")){ remove_tree($old) if -d $old; unlink($old) if -f $old; }
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
	my $swarmtar = "http://lotus2.earlham.ac.uk/lotus/packs/swarm2.1.13.zip";#"https://github.com/torognes/swarm/archive/master.zip";#"http://lotus2.earlham.ac.uk/lotus/packs/swarm206d.tgz";
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
	#vsearch
	print "Downloading vsearch executables..\n";
	my $vtars = "$bdir/vsearch.tar.gz";
	if ($isMac){
		#getS2("http://lotus2.earlham.ac.uk/lotus/packs/vsearch/vsearch-2.0.4-osx-x86_64/bin/vsearch",$vexe);
		getS2("https://github.com/torognes/vsearch/releases/download/v2.15.0/vsearch-2.15.0-macos-x86_64.tar.gz",$vtars);
	} else {
		#getS2("http://lotus2.earlham.ac.uk/lotus/packs/vsearch/vsearch-2.0.4-linux-x86_64/bin/vsearch",$vexe);
		getS2("https://github.com/torognes/vsearch/releases/download/v2.15.0/vsearch-2.15.0-linux-x86_64.tar.gz",$vtars);
	}
	run_cmd("tar", "-xzf", $vtars, "-C", $bdir);
	unlink($vtars) or warn "Could not remove $vtars: $!\n";
	my $vexe = $isMac ? "$bdir/vsearch-2.15.0-macos-x86_64/bin/vsearch" : "$bdir/vsearch-2.15.0-linux-x86_64/bin/vsearch";
	my ($vsearchVer,$vsearchStatus) = ("",-1);
	if (-s $vexe){
		run_cmd("chmod", "+x", $vexe);
		($vsearchVer,$vsearchStatus) = capture_cmd($vexe, "-v"); chomp $vsearchVer;
	}
	print "\n$vsearchVer\n";
	if (-s $vexe && $vsearchStatus == 0){# && $vsearchVer =~ m/vsearch v2.*/){ #not essential
		@txt = addInfoLtS("vsearch",$vexe,\@txt,1);
	} else {
		#system "rm $vexe";
		print "\n\nWARNING::\nvsearch exe did not exist at $vexe\n Therefore vsearch was not installed (fallback to usearch).\n\n";
		$finalWarning .= "vsearch exe did not exist at $vexe\n Therefore vsearch was not installed (fallback to usearch).\n";
	}

	#infernal
	print "Downloading infernal executables..\n";
	my $iexe = "$bdir/inf112.tar.gz";
	if ($isMac){
		getS2("http://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-macosx-intel.tar.gz",$iexe);
	} else {
		getS2("http://lotus2.earlham.ac.uk/lotus/packs/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz",$iexe);
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
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/VXtractor/vxtractor.pl",$vxexe);
	@txt = addInfoLtS("vxtractor",$vxexe,\@txt,1);
	$vxexe = "$bdir/vxtr/HMM.zip";
	getS2("http://lotus2.earlham.ac.uk/lotus/packs/VXtractor/HMMs.zip",$vxexe);
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
	my $fastt = "http://lotus2.earlham.ac.uk/lotus/packs/FastTree.c"; #http://www.microbesonline.org/fasttree/
	getS2($fastt,$exe1);
	$callret = system("gcc", "-DOPENMP", "-fopenmp", "-O3", "-finline-functions", "-funroll-loops", "-Wall", "-o", $exe, $exe1, "-lm");
	if ($callret != 0){
		print "\n\n=================\nProblem while compiling fasttree, trying fasttree without multithread and SSE support (might be slower, but if it's working..)\n";
		$finalWarning .= "fasttree compiled without multithreading support (you can not use the -thr LotuS option.\n";
		$exe = "$bdir/FastTree";
		$callret = system("gcc", "-DNO_SSE", "-O3", "-finline-functions", "-funroll-loops", "-Wall", "-o", $exe, $exe1, "-lm");}
	if ($callret != 0){
		$finalWarning .= "fasttree compilation failed. This is most likely an issue with your gcc version or the openMP libraries. See info on:\nhttp://www.microbesonline.org/fasttree/#Install\n";
		print "\n\n=================\nfasttree compilation failed. This is most likely an issue with your gcc version or the openMP libraries. See info on:\nhttp://www.microbesonline.org/fasttree/#Install\n"; exit(4);
	}

	run_cmd("chmod", "+x", $exe);
	@txt = addInfoLtS("fasttree",$exe,\@txt,1);


	#flash
	if (0){#not needed in lotus2 any longer..
		my $flashdir = $bdir."FLASH-1.2.10";
		my $fexe = "$flashdir/flash";
		my $tar = "$bdir/Flash.tar.gz";
		my $flashTar = "http://lotus2.earlham.ac.uk/lotus/packs/FLASH-1.2.10.tar.gz";#"http://sourceforge.net/projects/flashpage/files/FLASH-1.2.10.tar.gz/download";
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
	my $cdhitTar = "http://lotus2.earlham.ac.uk/lotus/packs/cd-hit_git.zip";#"https://github.com/weizhongli/cdhit/archive/master.zip";
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


	my $rdpf = "http://lotus2.earlham.ac.uk/lotus/packs/rdp_classifier_2.12.zip"; #"http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.6.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Frdp-classifier%2F&ts=1391590725&use_mirror=netcologne";
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
	my $clo = "http://lotus2.earlham.ac.uk/lotus/packs/clustalo-1.2.0-Ubuntu-x86_64";#"http://www.clustal.org/omega/clustalo-1.2.0-Ubuntu-x86_64";
	if ($isMac){
		$clo = "http://www.clustal.org/omega/clustal-omega-1.2.0-macosx";
	}
	$exe = $isMac ? "$bdir/clustal-omega-1.2.0-macosx" : "$bdir/clustalo-1.2.0-Ubuntu-x86_64";
	#system("wget -O $exe $clo");
	getS2($clo,$exe);
	run_cmd("chmod", "+x", $exe);
	@txt = addInfoLtS("clustalo",$exe,\@txt,1);
}

sub user_options(){

	if ($condaDBinstall){#no user input at all wanted
		return;
	}
	if ( $UID ne "??" || $forceUpdate || $usearchInstall ne ""){#a configured UID is sufficient to identify a previous installation
		my $inp="";
		
		if (!$forceUpdate && $usearchInstall eq ""){
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
		if ((0 && $inp eq "1") || $forceUpdate){ #normal online updater remains disabled; -forceUpdate is explicit
			my ($lsv,$msgEnd) = checkLtsVer($lver);
			#higher version? reinstall lotus3, autoinstall.pl, sdm
			if (version_is_newer($lsv,$lver) || $forceUpdate){
				print "New LotuS version available: updating from $lver to $lsv\n";
				my $updateArchive = "$ldir/files.tar.gz";
				my $updateDir = "$ldir/updates";
				my $installedHelper = "$ldir/helpers/autoInstall.pl";
				getS2("http://lotus2.earlham.ac.uk/lotus/lotus/updates/$lsv/files.tar.gz",$updateArchive);
				run_cmd("tar", "-xzf", $updateArchive, "-C", $ldir);
				die "Update archive is missing updates/autoInstall.pl\n" unless (-s "$updateDir/autoInstall.pl");
				die "Update archive is missing updates/lotus3\n" unless (-s "$updateDir/lotus3");
				if (-s $installedHelper != -s "$updateDir/autoInstall.pl" && !$forceUpdate){#at this point call autoupdate again
					print "Updated autoInstall.pl..\nAttempting to rerun autoInstall.pl\n";
					copy_file_atomic("$updateDir/autoInstall.pl", $installedHelper);
					exec($^X, $installedHelper, "-forceUpdate");
					die "Failed to rerun updated autoInstall.pl: $!\n";
				}
				for my $sourceDir (qw(sdm_src LCA_src rtk_src)){
					die "Update archive is missing updates/$sourceDir\n" unless (-d "$updateDir/$sourceDir");
				}
				# Compile the staged sources before replacing the installed source trees.
				my $nsdmp = compile_sdm("$updateDir/sdm_src");
				@txt = addInfoLtS("sdm",$nsdmp,\@txt,1);
				$nsdmp = compile_LCA("$updateDir/LCA_src");
				@txt = addInfoLtS("LCA",$nsdmp,\@txt,1);
				$nsdmp = compile_rtk("$updateDir/rtk_src");
				@txt = addInfoLtS("rtk",$nsdmp,\@txt,1);
				copy_file_atomic("$updateDir/autoInstall.pl", $installedHelper);
				copy_file_atomic("$updateDir/lotus3", "$ldir/lotus3");
				for my $sourceDir (qw(sdm_src LCA_src rtk_src)){
					replace_tree_atomic("$updateDir/$sourceDir", "$ldir/$sourceDir");
				}
				($lver,$sver) = getInstallVer("$ldir/sdm_src");
				remove_tree($updateDir) if -d $updateDir; unlink($updateArchive) if -e $updateArchive;
				if (length($msgEnd) >4){print "Additional information for this update:\n$msgEnd\n";}
				print "\nUpdated LotuS to version $lver\n\n";
				finishAI("u");
				exit(0);
			} else {
				print "You have the actual lotus version installed.\n"; exit(0);
			}
		} elsif($inp eq "2"){
			$onlyDbinstall = 1;
		}
	}
	#auto update END
	my $skipAll = 0 ; #debug option.. nerv
	if ($onlyDbinstall){
		print "Installing LotuS tax databases anew.. \nplease choose which databases to install in the following dialogs\n\n";
	}else{
		print "Total space required will be 0.3 - 5 GB.\nSome programs require a recent C++ compiler. Existing files are retained until their replacements download successfully, and lOTUs.cfg will be updated.\nContinue (y/n)?\nAnswer: ";
		my $confirmation = lc(read_user_input("installation confirmation"));
		if ($confirmation eq "y" || $confirmation eq "yes"){
			# continue
		} elsif ($confirmation eq "x") {
			$skipAll=1;
		} elsif ($confirmation eq "xx") {
			$skipAll=1;
			$refDBinstall[0] = 1;$refDBinstall[8] = 0;$ITSready=0;$getUTAX=0;
		} elsif ($confirmation eq "n" || $confirmation eq "no") {
			die "Installation cancelled by user.\n";
		} else {
			die "Invalid installation confirmation '$confirmation'.\n";
		}
		#print "\nThis is an experimental installer. Please send feedback and bug reports to: falk.hildebrand [at] gmail.com\n\n";
		if ($isMac){print "Mac system detected, installing corresponding mac software.\n";}

	#decide on blast
		if ($skipAll){
			return;
		}
		print "\n\nFor similarity based taxonomic assignments LotuS can either use \n (1) Blastn \n (2) Lambda \n (3) both, decide at runtime which to use or\n (0) none\n Answer:";
		while (1){
			my $choice = read_user_input("the similarity-search program choice");
			if ($choice =~ m/^[0123]$/){
				$installBlast = $choice;
				last;
			}
			print "Invalid answer; enter 0, 1, 2, or 3: ";
		}
	
	}


	#decide on database options

	print "\n\nDo you want to install a reference database 16S database for similarity based 16S annotations?\n";
	print " (1) KSGP (~1.5 GB), covering SSU for Archaea, Bacteria and Eukaryotes, 2026 release. \n (2) SILVA (~2.5 GB), contains LSU as well as SSU, 138.1 2020 release.\n (3) GreenGenes2 (~1 GB), 2022 release.\n (4) HITdb (~100 MB) 16S bacterial database specialized on the gut environment.\n";
	print " (5) PR2 (~100 MB), an SSU database specialized for marine eukaryotes.\n";
	print " (6) beeTax (~2 MB) database specialized (and named) on taxonomy specific to the bee gut.\n";
	print " (8) KSGP + SILVA + GG2 + PR2 + HITdb + beeTax (select a specific DB in each LotuS3 run)\n (0) no database.\n";
	print "Answer:";
	while (1){
		my $choice = read_user_input("the reference-database choice");
		if ($choice =~ m/^[0-6]$/ || $choice eq "8"){
			$refDBinstall[8] = 0;
			$refDBinstall[$choice] = 1;
			last;
		}
		die "Invalid reference-database choice '$choice'.\n";
	}
	#SILVA license
	if (!$skipAll && ($refDBinstall[2] || $refDBinstall[8])){
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

	print "\n\n -- ITS -- Do you want to\n (1) install databases and programs required to process ITS data (including fungi ITS UNITE database)\n (0) no ITS related packages\n Answer:";

	while (1){
		my $choice = read_user_input("the ITS package choice");
		if ($choice eq "1" || $choice eq "0"){
			$ITSready = $choice;
			last;
		}
		print "Invalid answer; enter 0 or 1: ";
	}

	#UTAX ref DBs..
	print "\n\n -- UTAX -- Do you want to\n (1) install utax taxonomic classification databases (16S, ITS)?\n (0) no utax related databases\n Answer:";
	while (1){
		my $choice = read_user_input("the UTAX database choice");
		if ($choice eq "1" || $choice eq "0"){
			$getUTAX = $choice;
			last;
		}
		print "Invalid answer; enter 0 or 1: ";
	}

}
