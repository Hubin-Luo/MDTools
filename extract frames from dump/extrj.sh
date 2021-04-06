#!/bin/bash
#
# Extract the trajectory of atoms in the simulation to an xyz file
# Execute the script followed by the traj file
#--------------------------------------------------------------
if [ $# -lt 1 ]; then
echo "input: trj_file [number of frame; default last]"
exit
fi
if [ -f extrj.xyz ]; then
rm extrj.xyz
fi

echo "#extracted from trj files" >> extrj.xyz
echo "    " >> extrj.xyz
if [ $# -gt 1 ]; then
  charline=`grep -n "ITEM: ATOMS" $1 | awk 'NR=='$2'{print}'`
else
  charline=`grep -n "ITEM: ATOMS" $1 | tail -1`
fi
nline=`echo $charline | awk -F: '{print $1}'`
id=`echo $charline | awk '{for(i=1;i<=NF;i++)
                              {
                               if($i=="id")
                                 {
                                  print i-2
                                  exit
                                 }
                              }
                           }'`
type=`echo $charline | awk '{for(i=1;i<=NF;i++)
                              {
                               if($i=="type")
                                 {
                                  print i-2
                                  exit
                                 }
                               }
                            }'`
x=`echo $charline | awk '{for(i=1;i<=NF;i++)
                               {
                                if($i=="x")
                                  {
                                   print i-2
                                   exit
                                  }
                                }
                          }'`
y=`echo $charline | awk '{for(i=1;i<=NF;i++)
                               {
                                if($i=="y")
                                  {
                                   print i-2
                                   exit
                                  }
                                }
                          }'`
z=`echo $charline | awk '{for(i=1;i<=NF;i++)
                               {
                                if($i=="z")
                                  {
                                   print i-2
                                   exit
                                  }
                                }
                          }'`
awk 'NR=='$nline'-5{print $1, "atoms"; print "XX atom types"; print "   "};
     NR=='$nline'-3{print $0, "xlo xhi"};
     NR=='$nline'-2{print $0, "ylo yhi"};
     NR=='$nline'-1{print $0, "zlo zhi"; print " "; print "Atoms"; print " "};
     NR >'$nline'{if($1=="ITEM:"){exit} else{print $'$id', $'$type', $'$x', $'$y', $'$z'}}' $1 >> extrj.xyz
nt=`cat extrj.xyz | awk 'BEGIN{n=1};
                           NR>=12{if($2>n)
                                   {
                                    n=$2
                                   }
                                 };
                           END{print n}'`
sed "s/XX/$nt/" extrj.xyz > temp
mv temp extrj.xyz
