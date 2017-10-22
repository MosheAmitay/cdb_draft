# ComplexDB

ComplexDB is a powerful tool, based on Python 3 and above, which enables searching heterodimers' complexes and antibodies' complexes with the RCSB Protein Data Bank (PDB). In addition, it's uniqueness is the ability to find the proteins that compound the chain's complex.
These functions use the PyPdB package.

The search is done according to specific criteria, some is determined by the user, and some-are default.
The finding is done by BLAST search. The user can filter the result, with specific criteria.

In addition to the name of the separate proteins, one can get the protein's length,the protein's sequence and the protein's identity to the chain of interest.

## Installation Requirements

In order to use the ComplexDB, install the PyPDB package:  [PyPDB](https://github.com/williamgilpin/pypdb)


## Usage

### Search

For searching, one must insert list of PDB entries, and search criteria, like number of chains in the complex and length of the chain. In antibody searching, one must supply also the maximal number of chain in the complex.

#### Complexes:

Get the complexes from a list of proteins (pdb codes), with given conditions and range

	Parameters:
	----------
	pdbs : list
	   A list where each value is a 4 character string giving a pdb entry of interest
	num_chains: int, #condition
	   The number of chains that the complex is constructed of.
	   Gets num_chains=2 as default
	length_chains: int, #condition
	   The minimum length of the chains that the complex is constructed of.
	   Gets num_chains=30 as default
	title_flag: bool, #condition
	   A boolean flag that define whether or not to check if the title of the PDB entry contains the words "Complex" or "Bound"
	   {Which is probably a complex)
	   Gets True as default
	counter: int, #range
	   The bottom edge of the range 
	   Gets counter=0 as default
	counter_end: int, #range
	   The top edge of the range 
	   Gets counter=0  as default in order to be replaced with len(pdbs) in the function
	Returns
	-------
	complexes : list
	   A list of potential string PDB 4 charachters entry according to the conditions that were given
	errors : list
	   A list of errors string PDB 4 charachters entry
	Examples
	--------
	s=Search()
	pdb=s.get_list_from_file()[0:500]
	s.search_complexes(pdbs)

#### Antibodies:

Get the antibodies from a list of proteins (pdb codes), with given conditions and range

	Parameters:
	----------
	pdbs : list
		A list where each value is a 4 character string giving a pdb entry of interest
	num_chains_min: int, #condition
		The minimum number of chains that the antibody is constructed of.
		Gets num_chains=2 as default
	num_chains_max: int, #condition
		The maximal number of chains that the antibody is constructed of.
		Gets num_chains=2 as default
	 length_chains: int, #condition
		   The minimum length of the chains that the complex is constructed of.
		   Gets num_chains=30 as default
	counter: int, #range
		The bottom edge of the range 
		Gets counter=0 as default
	counter_end: int, #range
		The top edge of the range 
		Gets counter=0  as default in order to be replaced with len(pdbs) in the function
	Returns
	-------
	antibodies : list
		A list of potential antibodies in string PDB 4 charachters entry according to the conditions that were given
	errors : list
		A list of errors string PDB 4 charachters entry
	Examples
	--------
	s=Search()
	pdb=s.get_list_from_file()[0:500]
	s.search_complexes(pdbs)
	
### Finding separate proteins:

For finding the separate proteins, one must supply the name of complex of interest (PDB code, 4 letters)

#### Complexes:

Find the separate proteins' names:


	Returns
	-------
	prot_names : list
		A list of PDB 4 characters entry, each represents seperate protein of chain.
   
	Examples
	--------
	c=Complex("1A0O")
	c.find_seperate_proteins()
	>>>['3CHY', '1FWP']
	
The details of the  separate proteins of a complex

	Returns
	-------
	self.__sep_prot_details : dictionary
	   A dictionary, where key is the seperate protein 4 charachters PDB name 
	   and value is tuple with the key's identity to the complex,sequence,flag-whether the seperate protein was found by accession number of the complex or not

	Examples
	--------
	c.get_proteins_details()
	>>>{'3CHY': (89, 'ADKELKFLVVDDFSTMRRIVRNLLKELGFNNVEEAEDGVDALNKLQAGGYGFVISDWNMPNMDGLELLKTIRADGAMSALPVLMVTAEAKKENIIAAAQAGASGYVVKPFTAATLEEKLNKIFEKLGM', False), '1FWP': (90, 'RQLALEAKGETPSAVTRLSVVAKSEPQDEQSRSQSPRRIILSRLKAGEVDLLEEELGHLTTLTDVVKGADSLSAILPGDIAEDDITAVLCFVIEADQITFETVEVSPKISTPPVLKLAAEQAPTGRVEREKTTR', False)}

Print the Complex, and it's details:
	
	Examples
	--------
	c=Complex("1A0O")
	print (c)
	>>>
	%%%%%%The Complex: 1A0O (Publication Year: 1998, Resolution: 2.95)%%%%%%
	%%% Protein: 3CHY %%%
	-Protein Name: CHEY
	-Protein Accession Number: P0AE67
	-Protein Chains: ACEG
	-Protein Chain Length: 128
	-Protein Sequence: ADKELKFLVVDDFSTMRRIVRNLLKELGFNNVEEAEDGVDALNKLQAGGYGFVISDWNMPNMDGLELLKTIRADGAMSALPVLMVTAEAKKENIIAAAQAGASGYVVKPFTAATLEEKLNKIFEKLGM
	-Protein Identity to Complex's Chain/s: 89
	-Protein Resolution: 1.66
	-Protein Publication Year: 1993
	%%% Protein: 1FWP %%%
	-Protein Name: CHEA
	-Protein Accession Number: P07363
	-Protein Chains: BDFH
	-Protein Chain Length: 134
	-Protein Sequence: RQLALEAKGETPSAVTRLSVVAKSEPQDEQSRSQSPRRIILSRLKAGEVDLLEEELGHLTTLTDVVKGADSLSAILPGDIAEDDITAVLCFVIEADQITFETVEVSPKISTPPVLKLAAEQAPTGRVEREKTTR
	-Protein Identity to Complex's Chain/s: 90
	-Protein Resolution: 0
	-Protein Publication Year: 1996
           

#### Antibodies:

Find the separate proteins' names:

	Returns
	-------
	prot_names : list
		A list of PDB 4 charachters entry, each represents seperate protein of chain.

	Examples
	--------
	b=Antibody("3TVM")
	b.find_seperate_proteins_anti()
	>>>['4ZAK', '3TVM']
	
The details of the  separate proteins of a complex

	Returns the details of the  seperate proteins of complex, refers to the complex

	Returns
	-------
	self.__sep_prot_details : dictionary
	   A dictionary, where key is the seperate protein 4 charachters PDB name 
	   and value is tuple with the key's identity to the complex,sequence,flag-whether the seperate protein was found by accession number of the complex or not

	Examples
	--------
	a=Antibody("1P2C")  
		  
	a.get_proteins_details()
	>>>{'4LZT': (100, 'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL', False), '2Q76': (89.0, 'EVQLEQSGAELMKPGASVKISCKATGYTFTTYWIEWIKQRPGHSLEWIGEILPGSDSTYYNEKVKGKVTFTADASSNTAYMQLSSLTSEDSAVYYCARGDGFYVYWGQGTTLTVSSASTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPWPSETVTCNVAHPASSTKVDKKIVPRDIELTQSPATLSVTPGDSVSLSCRASQSISNNLHWYQQKSHESPRLLIKYTSQSMSGIPSRFSGSGSGTDFTLSINSVETEDFGVYFCQQSGSWPRTFGGGTKLDIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRN', False)}
           
Print the Antibody, and it's details:

	a=Antibody("1P2C")
	print (a)
	>>>
	%%%%%%The Complex: 1P2C (Publication Year: 2004, Resolution: 2.00)%%%%%%
	%%% Protein: 4LZT %%%
	-Protein Name: Lysozyme C
	-Protein Accession Number: P00698
	-Protein Chains: CF
	-Protein Chain Length: 129
	-Protein Sequence: KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL
	-Protein Identity to Complex's Chain/s: 100
	-Protein Resolution: 0.95
	-Protein Publication Year: 1998
	%%% Antibody: 2Q76 %%%
	-Antibody Name: ,heavy chain VH+CH1 anti-lysozyme antibody F10.6.6,light chain anti-lysozyme antibody F10.6.6
	-Antibody Accession Number: ,P01868,P01837
	-Antibody Chains: BEAD
	-Antibody Chain Length: 430
	-Antibody Sequence: EVQLEQSGAELMKPGASVKISCKATGYTFTTYWIEWIKQRPGHSLEWIGEILPGSDSTYYNEKVKGKVTFTADASSNTAYMQLSSLTSEDSAVYYCARGDGFYVYWGQGTTLTVSSASTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPWPSETVTCNVAHPASSTKVDKKIVPRDIELTQSPATLSVTPGDSVSLSCRASQSISNNLHWYQQKSHESPRLLIKYTSQSMSGIPSRFSGSGSGTDFTLSINSVETEDFGVYFCQQSGSWPRTFGGGTKLDIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRN
	-Antibody Identity to Complex's Chain/s: 89.0
	-Antibody Resolution: 2.00
	-Antibody Publication Year: 2008
