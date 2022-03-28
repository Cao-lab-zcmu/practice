read -p "The Display Mode >>> " parameter
if [ $parameter == 0 ]
then
  i=0
elif [ $parameter == 1 ]
then
  i=9
fi
xrandr -s $i
