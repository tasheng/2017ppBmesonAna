i=0
f=1000




PtMin=7
PtMax=10

while [ $i -lt $f ]

do 
	echo "Now Working on i = " $i 

	source doRoofit.sh $i  ${PtMin} ${PtMax}

	i=$(($i+1))

done

