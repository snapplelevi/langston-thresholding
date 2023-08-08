all:
	echo "test"

.PHONY: clean
clean:
		rm ./src/*.o
		rm ./src/thresholding.dll