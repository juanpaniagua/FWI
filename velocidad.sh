#! /bin/sh
# File: model2.sh

# Set messages on
set -x



# Name output binary model file
modfile=velocity.dat

# Name output encapsulated Postscript image file
psfile=velocity.eps

# Remove previous .eps file
rm -f $psfile

trimodel xmin=0 xmax=16 zmin=0 zmax=16 \
1 xedge=0,16 \
  zedge=0,0 \
  sedge=0,0 \
2 xedge=0,16\
  zedge=3.0,10.0 \
  sedge=0,0 \
3 xedge=0,16 \
  zedge=16,16 \
  sedge=0,0 \
 kedge=1,2,3 \
 sfill=1,0.1,0,0,0.12,0,0 \
 sfill=1,4.5,0,0,0.08,0,0 > $modfile
##       x,z

# Create a Postscript file of the model
spsplot < $modfile > $psfile \
        gedge=0.5 gtri=3.0 gmin=0.2 gmax=0.8 \
titlecolor=black axescolor=black \
        labelz="Depth z (m)" labelx="Distance x (m)" \
        dxnum=1.0 dznum=1.0 box=16.0 hbox=8.0
        
        
echo "Generando modelo de velocidades triangulado."

tri2uni < velocity.dat n1=500 n2=500 d1=0.032 d2=0.032 > velocity.bin

# Exit politely from shell
exit

