{
cs=($5 * $6);
ts=($5 * $7);
if ($4 == "F") str= "+";
else if ($4 == "R") str= "-";

print $2"\t"$3"\t"$3+1"\t"int(cs)"\t"int(ts)"\t"str;
}
