# sequenceAlignment-Bioinformatics

Implementation (python 3.4) of global sequence alignment based on dynamic programming approach.  </br>

Usage:  </br>

python sequenceAlignment.py -i1 (seq1 file) -i2 (seq2 file) -s (scoring file) -g (gap penalty)  </br>

Score file example:  </br>

	A	C	G	T	</br>
A	91	-114	-31	-123   </br>
C	-114	100	-125	-31    </br>
G	-31	-125	100	-114   </br>
T	-123	-31	-114	91     </br>
