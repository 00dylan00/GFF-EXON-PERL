# GFF-EXON
Functions:
	(1.1) Hash of length of all transcripts.
		Open Gff File and create a Hash of Arrays. The first key of the hash
		corresponds to the ID of the gene the second key to the transcript ID
		and the values are an array of the transcript lengths.
	(1.2) Hash of total exon length.
		Transverse the hash and create a Hash of Arrays. The first key of the hash
		corresponds to the ID of the gene the second key to the transcript ID
		and the values are the sum total of elements in the array (total exon
		length).
	(1.3) Hash of all GFF features.
		Open GFF File and create a Hash of Arrays. For this hash the first key
		corresponds to the ID of the transcripts and the values of the array are
		the GFF lines corresponding to said tranacript.
	(1.4) Shortest transcript features.
		Traverse the Hash of total exon length and determine the shortest exon  
		transcript for each gene. Once found, traverse the Hash of all features  
		and retrieve all features corresponding to said transcript.
	(1.5) Longest transcript features.
		Traverse the Hash of total exon length and determine the longest exon  
		transcript for each gene. Once found, traverse the Hash of all features  
		and retrieve all features corresponding to said transcript.
	(1.6) Sort transcripts in descending order.
		Traverse the Hash of total exon length and sort the secondary keys in 
		descending order.
	(1.7) Sort transcripts in ascending order.
		Traverse the Hash of total exon length and sort the secondary keys in 
		ascending order.  
	(1.8) Shortest Gene transcript.
 		Traverse the Hash of total exon length and find the shortest exon in 	
		the whole GFF file. Output gene and transcript ID. Option to print 
		all GFF features associated to transcript depending command line 
		response. Traverse Hash of all GFF features to retrieve information.
	(1.9) Longest Gene transcript.
 		Traverse the Hash of total exon length and find the longest exon in 	
		the whole GFF file. Output gene and transcript ID. Option to print 
		all GFF features associated to transcript depending command line 
		response. Traverse Hash of all GFF features to retrieve information.
	(1.10) Statistical analysis.
		Carry out statistical analysis of all exons in GFF file to calculate
		number of total exons, mean, standard deviation and variance. Traverse
		Hash of total exons and push secondary key values into an array. 
	(1.11) Plot differences of longest and shortest exon for each gene.
		Traverse hash of exons and calculate shortest and longest exon for
		each. Compute difference and represent value in logarithmic scale
		base 2. Only output genes with different exon length. 
	(1.12) Search for gene.
		Search for gene of interest by traversing Hash of exons and outputting
		exon length. Option to print all GFF features associated to transcript 
		depending command line response. Traverse Hash of all GFF features to 
		retrieve information. This works because transcript ID starts with 
		gene ID. 
