#!/usr/bin/perl/ -w
$|++;
use strict;
use File::Path;
use Time::HiRes qw( time sleep );

######################################################################################################################################################
#
#	Description
#		This is a perl script to read the jnctnInfo.txt generated from HMMSplicerBEDToSAMParser.pl, and select the jucntions based on user defined score
#	and total supported lengths cut-offs. It will run weblog for the jucntion flanking sequence, and determine if there are conserved site. 
#
#	Input	
#
#		--readCutoff=	the minimum number of support reads for junctions, default = 10;
#		--scoreCutoff=	the minimum score for junctions to be analyzed, default = 600;
#		--jnctnInfoPath= the uncollapsed BED file contains the BED info for all unfiltered reads;
#		--entropyCutoff= the minimum ratio of entropy to that of the spliced sites to be considered as conserved; default=0.5;
#		--minConNum=	the minimum number of consevered residues;
#		--trimExonEnd=	the number of bases to be trimmed from the exon end of the flanking seq, default = 20;
#		--trimIntronEnd= the number of bases to be trimmed from the intron end of the flanking seq, default = 20;
#		--minJunctNum=		junctionType with less than this number of junctions will not be plotted;
#		--minJunctDis=		junctions with 5' donor site closer than this distance is regarded as potential artifact of internal repeats and will not be counted in the same se junct type;
#		--outDir= the directory for output;
#
#	Output
#
#	Usage
#		perl junctSeqConservationAnalyzer_v0.1.pl --jnctnInfoPath=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/HMMSplicerBEDToSAMParser/HMMSplicer/junctionInfo/HMMSplicer.filter.noncan.unique.JnctnInfo.txt --trimExonEnd=20 --trimIntronEnd=20 --readCutoff=5 --scoreCutoff=900 --entropyCutoff=0.5 --minConNum=2 --minJunctNum=20 --minJunctDis=200 
#
#	Assumption 
#
#	Versions
#
#		v0.1
#		-debut;
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#


#1----------Read the parameters----------#
use vars qw ($readCutoff $scoreCutoff $jnctnInfoPath $entropyCutoff $minConNum $trimExonEnd $trimIntronEnd $minJunctNum $minJunctDis $outDir $paraTag);
my ($readCutoff, $scoreCutoff, $jnctnInfoPath, $entropyCutoff, $minConNum, $trimExonEnd, $trimIntronEnd, $minJunctNum, $minJunctDis, $outDir, $paraTag) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

my ($flankSeqByJunctSeqHsh_ref, $allInfoHash_ref) = readJunctionInfo();

calculateEntropyAndPlotLogo($flankSeqByJunctSeqHsh_ref, $allInfoHash_ref);
checkConservedResidue();

printCMDLogOrFinishMessage("finishMessage");
exit;
#========================================================= Main body ends ===========================================================================#
################################################################## readParameters ####################################################################
sub readParameters {
	
	$readCutoff = 10;
	$scoreCutoff = 1000;
	$entropyCutoff = 0.5;
	$minConNum = 2;
	$trimExonEnd = 20;
	$trimIntronEnd = 20;
	$minJunctNum = 20;
	$minJunctDis = 200;
	$outDir = "./";

	foreach (@ARGV) {
		if ($_ =~ m/--readCutoff=/) {$readCutoff = substr ($_, index ($_, "=")+1);} 
		elsif ($_ =~ m/--scoreCutoff=/) {$scoreCutoff = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--jnctnInfoPath=/) {$jnctnInfoPath = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--entropyCutoff=/) {$entropyCutoff = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--minConNum=/) {$minConNum = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--trimExonEnd=/) {$trimExonEnd = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--trimIntronEnd=/) {$trimIntronEnd = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--minJunctNum=/) {$minJunctNum = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--minJunctDis=/) {$minJunctDis = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--outDir=/) {$outDir = substr ($_, index ($_, "=")+1);}
	}
		
	my $paraTag = "r$readCutoff.s$scoreCutoff.e$entropyCutoff.c$minConNum.j$minJunctNum.m$minJunctDis";
	$outDir =~ s/\/$// if ($outDir =~ m/\/$/);
	system ("mkdir -p -m 777 $outDir/$paraTag/weblogo/"); #---to make sure the output
	system ("mkdir -p -m 777 $outDir/$paraTag/entropy/"); #---to make sure the output
	system ("mkdir -p -m 777 $outDir/$paraTag/seq/"); #---to make sure the output
	
	open (TEST, $jnctnInfoPath) || die "Cannot open $jnctnInfoPath\n"; close TEST;
	
	return ($readCutoff, $scoreCutoff, $jnctnInfoPath, $entropyCutoff, $minConNum, $trimExonEnd, $trimIntronEnd, $minJunctNum, $minJunctDis, $outDir, $paraTag);
}
################################################################## readJunctionInfo ####################################################################
sub calculateEntropyAndPlotLogo {
	
	my %flankSeqByJunctSeqHsh = %{$_[0]};
	my %allInfoHash = %{$_[1]};
	
	my $junctTypeNum = keys %flankSeqByJunctSeqHsh;
	
	printProgressScale("Calculating the entropy of the aligned positions and generating weblogo", 50);
	my $progCount=0;
	my $allPdfPath = "";
	my $errorLogPath = "$outDir/$paraTag/error.log.txt";
	open (ERROR, ">$errorLogPath");
	close ERROR;
	my $filterSeqLogPath = "$outDir/$paraTag/filtered.JunctInfo.txt";
	open (FILTERINFO, ">$filterSeqLogPath");
	foreach my $junctSeqType (sort {$b cmp $a} keys %flankSeqByJunctSeqHsh) {
		$progCount++;

		my $entropyPath = "$outDir/$paraTag/entropy/$junctSeqType.entropy.txt";
		my $seqPath = "$outDir/$paraTag/seq/$junctSeqType.seq.fasta";
		my $logoPath = "$outDir/$paraTag/weblogo/$junctSeqType.weblogo.pdf";
		
		$allPdfPath .= " \'$logoPath\'";

		my $length;
		
		open (SEQFILE, ">$seqPath");
		my %strndCountHsh;
		foreach my $junctStr (sort {$a cmp $b} keys %{$flankSeqByJunctSeqHsh{$junctSeqType}}) {
			my $originalInfoLine = $allInfoHash{$junctStr};
			print FILTERINFO $originalInfoLine."\n";
			my $seq = ${$flankSeqByJunctSeqHsh{$junctSeqType}}{$junctStr};
			$length = length $seq;
			print SEQFILE ">".$junctStr."\n";
			print SEQFILE $seq."\n";
			my $strnd = chop $junctStr;
			$strndCountHsh{$strnd}++;
		}
		close SEQFILE;
		my $seqNum = keys %{$flankSeqByJunctSeqHsh{$junctSeqType}};
		
		my $title = $junctSeqType.".".$seqNum;
		updateProgressBar("$junctSeqType", $progCount, $junctTypeNum, 50, 5);
		my $weblogoEntropyCMD = "weblogo -f $seqPath -F txt -n $length -s large --title $title -A rna -c classic >$entropyPath";
		my $grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
		runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outDir/error.log.txt");

		$weblogoEntropyCMD = "weblogo -f $seqPath -F pdf -n $length -s large --title $title -A rna -c classic >$logoPath";
		$grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
		runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, $errorLogPath);

	}
	
	print "\n\n";
	close FILTERINFO;
	close SEQFILE;
	
	print "Merging all pdfs.\n";
	my $mergePdfCMD = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=\'$outDir/$paraTag/all.merged.pdf\' $allPdfPath";
	my $grepCMD = "ps -ef | grep $allPdfPath | grep -v grep";
	runAndCheckSerialTask($grepCMD, $mergePdfCMD, $mergePdfCMD, $errorLogPath);
	
}
############################################################ IGVCovSub #####################################################################
sub runAndCheckSerialTask {

	my $grepCmd = $_[0];
	my $grepStr = $_[1];
	my $cmd = $_[2];
	my $errorLogPath = $_[3];

	system (qq|$cmd 2>>$errorLogPath &|);
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCmd`;
		sleep (0.001);
	}
}
#####################################updateProgressBar###########################################################################################
sub updateProgressBar {
	
	my $strToPrint = $_[0];
	my $progressCount = $_[1];
	my $totalCount = $_[2];
	my $scaleMax = $_[3];
	my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
	my $progressPct = int (($progressCount/$totalCount)*$scaleMax);

	my $progressBar = "|";
	for my $i (1..$progressPct) {$progressBar .= ">";}
	for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
	$progressBar .= "|";

	my $extraWhiteSpaceStr = "";
	for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
	print $progressBar.$strToPrint.$extraWhiteSpaceStr."\r";

}
#####################################print Progress Scale###########################################################################################
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
################################################################## readJunctionInfo ####################################################################
sub readJunctionInfo {
				
	#non-cannoical	unique	DS571145,107720,107770	1276.7	13	5	8	TCAAGGAAATACAGGTTCAACAACTGAAGGTTTGTTTCATTAAATAGTGTATAAATAATTCA	TATAAATAATTCACGATATTCATCTTTTAAGAATGTACCGTTTCATGGGGAACTTGTAATCA	TTGA	*	unq17	CTTCAAGTAGACACATGTGCTATTTGTAGAAATAGTCTTATGGAATTATGTTTAGAATGTCAAGGAAATACAGGTTCAACAACTGAAGG	ATGTACCGTTTCATGGGGAACTTGTAATCATGCATTCCATACTCATTGTATTTCTTCATGGTTAAGACAACGTGCTGTCT
	#non-cannoical	unique	DS571145,107722,107772	852.4	1	0	1	AAGGAAATACAGGTTCAACAACTGAAGGTTTGTTTCATTAAATAGTGTATAAATAATTCACG	TAAATAATTCACGATATTCATCTTTTAAGAATGTACCGTTTCATGGGGAACTTGTAATCATG	TGAT	*	unq18	AAATAGTCTTATGGAATTATGTTTAGAATGTCAAGGAAATACAGGTTCAACAACTGAAGGTT	GTACCGTTTCATGGGGAACTTGTAATCATGCATTCCAT
	#non-cannoical	unique	DS571145,119991,120040	858.0	2	1	1	ATTTCTTTGATCAAATCATCTTCTTGTTTTTTCTATTAAAATAATTAATAATTTTTCTAAAA	TAATTTTTCTAAAAAAAATATATGAAACAAACTATTTCTTGATATTCCTTCATCAAACCTTC	TTAC	*	unq30	CAATTTCTTTGATCAAATCATCTTCTTGTTTT	TATTTCTTGATATTCCTTCATCAAACCTTCGAATGATTGTTCTGTTTTTGTTTCATTTTTACTCTCATTTCTATCTTGC
	#non-cannoical	unique	DS571145,119992,120038	843.4	1	0	1	TTTCTTTGATCAAATCATCTTCTTGTTTTTTCTATTAAAATAATTAATAATTTTTCTAAAAA	AATAATTTTTCTAAAAAAAATATATGAAACAAACTATTTCTTGATATTCCTTCATCAAACCT	TCAA	*	unq32	GTATTTACATCATCAACTCCTTTGACCATCTCTTCAATTTCTTTGATCAAATCATCTTCTTGTTTTT	ACTATTTCTTGATATTCCTTCATCAAACCTTCG
	#non-cannoical	unique	DS571145,119992,120040	1070.5	2	1	1	TTTCTTTGATCAAATCATCTTCTTGTTTTTTCTATTAAAATAATTAATAATTTTTCTAAAAA	TAATTTTTCTAAAAAAAATATATGAAACAAACTATTTCTTGATATTCCTTCATCAAACCTTC	TCAC	*	unq33	CCATCTCTTCAATTTCTTTGATCAAATCATCTTCTTGTTTTT	TATTTCTTGATATTCCTTCATCAAACCTTCGAATGATTGTTCTGTTTTTGTTTCATTTTTACTCTCATTTCTA

	my (%flankSeqByJunctSeqHsh, %allInfoHash);
	my $validJunctNum = 0;
	
	open (JUNCTINFO, $jnctnInfoPath);
	while (my $theLine = <JUNCTINFO>) {
		chomp $theLine;
		my @theLineSplt = split /\t/, $theLine;
		my $score = $theLineSplt[3];
		my $supportRead = $theLineSplt[4];
		
		if (($score >= $scoreCutoff) and ($supportRead >= $readCutoff)) {
			
			my $leftFlankSeq = $theLineSplt[7];
			my $rightFlankSeq = $theLineSplt[8];
			my $junctStr = $theLineSplt[2];
			my $oriLen = length $leftFlankSeq;
			my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
			die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0);
			
			my $trimLeftFlankSeq = substr $leftFlankSeq, $trimExonEnd, $trimLen;
			my $trimRightFlankSeq = substr $rightFlankSeq, $trimIntronEnd, $trimLen;
			my $concateSeq = $trimLeftFlankSeq.$trimRightFlankSeq;
			
			my $junctSeq = my $revComJunctSeq = $theLineSplt[9];
			$revComJunctSeq = reverse $revComJunctSeq;
			$revComJunctSeq =~ tr/ACGTacgt/TGCAtgca/;
			my $strnd = "+";

			#---check if the reverse complement is present, if yes, reverse everything 
			if (exists $flankSeqByJunctSeqHsh{$revComJunctSeq}) {
				$junctSeq = $revComJunctSeq;
				$strnd = "-";
				$concateSeq = reverse $concateSeq;
				$concateSeq =~ tr/ACGTacgt/TGCAtgca/;
			}

			#---convert to RNA
			$concateSeq =~ tr/Tt/Uu/;
						
			$junctStr = $junctStr.",".$strnd;
			
			$allInfoHash{$junctStr} = $theLine;

			if (not exists $flankSeqByJunctSeqHsh{$junctSeq}) {
				
				${$flankSeqByJunctSeqHsh{$junctSeq}}{$junctStr}=$concateSeq;

			} else {
			
				my $withinMinDis = "no";
				my @junctStrSplt = split /\,/, $junctStr;
				my $cntg = $junctStrSplt[0];
				my $junctMidPt = int ($junctStrSplt[1] + ($junctStrSplt[2]-$junctStrSplt[1])/2);

				foreach my $storedJunctStr (keys %{$flankSeqByJunctSeqHsh{$junctSeq}}) {
					my @storedJunctStrSplt = split /\,/, $storedJunctStr;
					my $storedCntg = $storedJunctStrSplt[0];
					my $storedJunctMidPt = int ($storedJunctStrSplt[1] + ($storedJunctStrSplt[2]-$storedJunctStrSplt[1])/2);
					if ($cntg eq $storedCntg) {#--on the same contig and 
						my $midPtDis = abs ($storedJunctMidPt - $junctMidPt);
						$withinMinDis = "yes" if ($midPtDis < $minJunctDis);
					}
				}
				
				${$flankSeqByJunctSeqHsh{$junctSeq}}{$junctStr}=$concateSeq if ($withinMinDis eq "no");
			}

			my $junctSeqNum = keys %flankSeqByJunctSeqHsh;
			$validJunctNum++;
			
			print "$junctSeqNum types of junctions and $validJunctNum junctions were stored.\r";
		}
	}
	close JUNCTINFO;
	
	print "\n\n";
	
	#---remove the junctSeqType with < 3 seq
	my $rmNum = 0;
	foreach my $junctSeqType (keys %flankSeqByJunctSeqHsh) {
		my $seqNum = keys %{$flankSeqByJunctSeqHsh{$junctSeqType}};
		if ($seqNum < $minJunctNum) {
			delete $flankSeqByJunctSeqHsh{$junctSeqType};
			$rmNum++;
			print "$junctSeqType type junctions has less than $minJunctNum junctions. $rmNum types Removed.\r";
		}
	}
	print "\n\n";

	return (\%flankSeqByJunctSeqHsh, \%allInfoHash);

}
#################################################### print command log ###############################################################################
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
#################################################### print command log ###############################################################################
sub checkConservedResidue {
	
	my $allEntropyFileHsh_ref = getFilePathInDir("$outDir/$paraTag/entropy/", "\.entropy\.txt");
	my %allEntropyFileHsh = %{$allEntropyFileHsh_ref};
	my %cnsvJunctTypeHsh;
	
	foreach my $fileName (sort {$a cmp $b} keys %{$allEntropyFileHsh_ref}) {
		my $filePath = $allEntropyFileHsh{$fileName};
		my @fileNameSplt = split /\./, $fileName;
		my $junctSeqType = $fileNameSplt[0];
		open (ENTROPYFILE, "$filePath") || die "Cannot open $filePath :$!";
		my %posEntropyHsh;
		my $maxEntropy = 0;
		while (my $theLine = <ENTROPYFILE>) {
			chomp $theLine;
			next if (($theLine =~ m/^\#/) or ((length $theLine) < 1));
			my @theLineSplt = split /\t/, $theLine;
			my $pos = $theLineSplt[0];
			my $entropy = $theLineSplt[5];
			$posEntropyHsh{$pos} = $entropy;
			$maxEntropy = $entropy if ($entropy > $maxEntropy);
		}
		close ENTROPYFILE;
		
		my $cnsvResidue = 0; #----count the 4 100% conserved junctSeq
		my $ratioCutoff = $maxEntropy*$entropyCutoff;

		foreach my $pos (keys %posEntropyHsh) {
			$cnsvResidue++ if ($posEntropyHsh{$pos} >= $ratioCutoff);
		}
		$cnsvResidue = $cnsvResidue - 4;
		$cnsvJunctTypeHsh{$junctSeqType} = $cnsvResidue;
	}
	
	my $cnsvResidueLogPath = "$outDir/$paraTag/conserved.residueNum.log.txt";
	my $fileToMerge = "";
	open (CNSVRESIDUE, ">$cnsvResidueLogPath");
	foreach my $junctSeqType (sort {$cnsvJunctTypeHsh{$b} <=> $cnsvJunctTypeHsh{$a}} keys %cnsvJunctTypeHsh) {
		my $filePath = "$outDir/$paraTag/weblogo/$junctSeqType.weblogo.pdf";
		print CNSVRESIDUE $junctSeqType."\t".$cnsvJunctTypeHsh{$junctSeqType}."\n";
		if ($cnsvJunctTypeHsh{$junctSeqType} >= $minConNum) {
			$fileToMerge .= " $filePath";
		}
		
	}
	close CNSVRESIDUE;
	
	print "Merging the pdf of conserved alignements.\n";
	my $errorLogPath = "$outDir/$paraTag/error.log.txt";
	my $mergePdfCMD = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=\'$outDir/$paraTag/conserved.merged.pdf\' $fileToMerge";
	my $grepCMD = "ps -ef | grep $fileToMerge | grep -v grep";
	runAndCheckSerialTask($grepCMD, $mergePdfCMD, $mergePdfCMD, $errorLogPath);
}
#################################################### print command log ###############################################################################
sub getFilePathInDir{

	my $pathToRead = $_[0];
	my $strToCheck = $_[1];
	
	my %filePathHsh;
	#---add the back slash if the path does not contain one
	$pathToRead = $pathToRead."/" if ($pathToRead !~ m/\/$/);
	opendir(DIRHANDLE, "$pathToRead") || die "Cannot read $pathToRead: $!";
	foreach my $fileName (sort {$b cmp $a} readdir(DIRHANDLE)) {
		if ($fileName =~ m/$strToCheck/){
			my $fullPath = $pathToRead.$fileName;
			$filePathHsh{$fileName} = $fullPath;
		}
	}
	return \%filePathHsh;
	close DIRHANDLE;
}
