#!/usr/bin/perl

open(FP,"$ARGV[0]");
$var="$ARGV[0]";
#print $ARGV[0];
#print $var;
$var=~s/\.temp//g;
#print $var;
print"\@relation $var\n\n";
while($line5=<FP>)
{
if($line5=~/ID/)
{
@q=split"\t",$line5;
for($o=2;$o<@q;$o++)
{
	chomp $q[$o];
print "\@attribute $q[$o] numeric\n";
}
print"\@attribute class numeric\n";
print"\@data\n";
}
else
{
#$line=~s/\d+.mol2\t//g;
@w=split"\t",$line5;
for($o=2;$o<@w;$o++)
{
chomp($w[$o]);
print"$w[$o],"; 
}
chomp($w[1]);
print"$w[1]\n";
}
}
