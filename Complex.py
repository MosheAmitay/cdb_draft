from pypdb.pypdb import *
from IComplex import IComplex
from Search import *
class Complex(IComplex):

    '''
    =================
    Functions for for the user
    =================
    '''

    def __init__(self,name_pdb):
        '''Initiating a Complex object, with given PDB code. Using the IComplex initiation.
           Parameters
           ----------
           name_pdb : string
               A 4 character string giving a pdb entry of interest
           
           '''
        IComplex.__init__(self, name_pdb)
        self.__sep_prot_details={}

    def find_seperate_proteins(self):
        '''Find the seperate proteins of complex, using find_protein function
    
        Returns
        -------
        prot_names : list
            A list of PDB 4 charachters entry, each represents seperate protein of chain.
       
        Examples
        --------
        c=Complex("1A0O")
        c.find_seperate_proteins()
        >>>['3CHY', '1FWP']

        '''
        prot_names=[]
        for chain_index in range(len(self.all_info["polymer"])):
            pdb_code,identity,seq,flag=self.find_protein(chain_index=chain_index)
            if pdb_code=="None":
                self.__sep_prot_details={}
                self.__sep_prot_details["None"]="none"
                return []
            this_details=(identity,seq,flag)
            self.__sep_prot_details[pdb_code]=this_details
            prot_names.append(pdb_code)
        return prot_names

    def get_proteins_details(self):
        '''Returns the details of the  seperate proteins of complex, refers to the complex

           Returns
           -------
           self.__sep_prot_details : dictionary
               A dictionary, where key is the seperate protein 4 charachters PDB name 
               and value is tuple with the key's identity to the complex,sequence,flag-whether the seperate protein was found by accession number of the complex or not

           Examples
           --------
           c.get_proteins_details()
            >>>{'3CHY': (89, 'ADKELKFLVVDDFSTMRRIVRNLLKELGFNNVEEAEDGVDALNKLQAGGYGFVISDWNMPNMDGLELLKTIRADGAMSALPVLMVTAEAKKENIIAAAQAGASGYVVKPFTAATLEEKLNKIFEKLGM', False), '1FWP': (90, 'RQLALEAKGETPSAVTRLSVVAKSEPQDEQSRSQSPRRIILSRLKAGEVDLLEEELGHLTTLTDVVKGADSLSAILPGDIAEDDITAVLCFVIEADQITFETVEVSPKISTPPVLKLAAEQAPTGRVEREKTTR', False)}

           '''
        if len(self.__sep_prot_details.keys()) == 0:  # the separated proteins wasn't found yet
            not_to_use = self.find_seperate_proteins()

        return self.__sep_prot_details

    def get_complex_details(self, index, protein):
        '''Get the complexes deatails

          Parameters:
          ----------
          index : int
              The index of the chain of indetest
          protein: string
              A 4 charachters string, represents the PDB entry for the seperate protein
          Returns
          -------
          details : list
              A list wit protein's details: name, accession numer,chains were checked, chains length,resolution,publication year
     
          Examples
          --------
          c=Complex("1A0O")
          name,identity,seq,bal=c.find_protein(chain_index=1)
          c.get_complex_details(1,name)
          >>>['CHEA', 'P07363', 'BDFH', 134, 0, '1996']

          '''

        details = []
        description = describe_pdb(protein)

        chain_id, chain_length, acc_num = self.get_search_params(index)
        name = self.all_info["polymer"][index]["polymerDescription"]["@description"]
        resolution = self.set_resolution(description)
        pub_year = description['release_date'][0:4]
        details.append(name)
        details.append(acc_num)
        details.append(chain_id)
        details.append(chain_length)
        details.append(resolution)
        details.append(pub_year)
        return details

    def __str__(self):
        '''Print function
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
           
          '''
        if len(self.__sep_prot_details.keys())==0: #the separated proteins wasn't found yet
            not_to_use=self.find_seperate_proteins()
        if "None" in self.__sep_prot_details.keys():
            return "There is no seperate protein for this complex"
        details_print="%%%%%%"+IComplex.__str__(self)+"%%%%%%\n"

        index=0
        for key,value in self.__sep_prot_details.items():
            newkey=key
            details = self.get_complex_details( index,key)
            if value[2]==True: newkey+="*"
            details_print+="\n%%% Protein: "+newkey+" %%%"
            details_print+="\n-Protein Name: "+details[0]
            details_print+="\n-Protein Accession Number: "+details[1]
            details_print += "\n-Protein Chains: "+details[2]
            details_print += "\n-Protein Chain Length: "+str(details[3])
            details_print += "\n-Protein Sequence: "+value[1]
            details_print += "\n-Protein Identity to Complex's Chain/s: "+str(value[0])
            details_print += "\n-Protein Resolution: "+str(details[4])
            details_print += "\n-Protein Publication Year: "+details[5]+"\n"
            index+=1
        return details_print



    def find_protein(self,chain_index=0,identity_percent=80,limit_percent=20):
        '''Find the seperate protein of a given chain index, with given conditions

          Parameters:
           ----------
           chain_index : int
               The index of chain of interest, to get chain's details
               Gets chain_index=0 as default
           identity_percent: int, #condition
               The minimal percent of identity to filter.
               Gets identity_percent=80 as default
           limit_percent: int, #condition
               The maximal percent of chain length, to be considered as identical to the given chain
               Gets limit_percent=29 as default
           Returns
           -------
           pdb_code : string
                A 4 character string represent a PDB entry of the seperate protein
           identity : int
               The identity percent of the seperate protein, calculated with BLAST
            seq : string
                The sequence of the chain of interest
            flag : bool
               Boolean variable, represents whether the seperate protein was found by accession number of the complex or not
           Examples
           --------
            c=Complex("1A0O")
            c.find_protein(chain_index=1)
            >>>('1FWP', 90, 'RQLALEAKGETPSAVTRLSVVAKSEPQDEQSRSQSPRRIILSRLKAGEVDLLEEELGHLTTLTDVVKGADSLSAILPGDIAEDDITAVLCFVIEADQITFETVEVSPKISTPPVLKLAAEQAPTGRVEREKTTR', False)

           '''
        chain_id, chain_length, acc_num = self.get_search_params(chain_index)

        blast_results = get_blast2(self.pdb_id, chain_id=chain_id)
        results = []

        # getting the data from the blast: pdb_name,eval,identity and length of chain
        for res in range(len(blast_results[1])):
            this_result = self.__get_hit_params(blast_results[1][res])
            if this_result[2] == "": #not seperate protein for sure
                continue
            if int(this_result[2]) >= identity_percent:
                if this_result[0] == self.pdb_id: continue
                all_info = get_all_info(this_result[0])
                if type(all_info["polymer"]) is dict: #check this is seperate protein
                    tup1 = (this_result[0], int(this_result[2]), int(this_result[1]), this_result[3])
                    results.append(tup1)

        # if there is no identity>80%:
        if len(results) == 0:
            return "None", "None", "None", False

        # sort the list according to length
        sorted_result = sorted(results, key=lambda x: x[2], reverse=True)

        # the identity is limited to 80, therefore the limit of length is 20% of the given length
        limit = int(chain_length / 100) * limit_percent
        len_range = 0
        i = 0
        length_result = []
        flag = False
        # find the closest length to the chain_length
        while len(length_result) == 0 and len_range <= limit:
            i = 0
            while i < len(sorted_result) and sorted_result[i][2] > chain_length + len_range:
                i += 1
            # sorted_result[i][2]>=chain_length+len_range and
            while i < len(sorted_result) and sorted_result[i][2] <= chain_length + len_range and sorted_result[i][
                2] >= chain_length - len_range:
                length_result.append(sorted_result[i])
                i += 1
                if i == len(sorted_result) or not (sorted_result[i][2] <= chain_length + len_range and sorted_result[i][
                    2] >= chain_length - len_range):
                    flag = True

            if flag == True:
                break
            len_range += 1

        # no protein in the limited length
        if len(length_result) == 0:
            if acc_num=="":
                return "None", "None", "None", False
            for k in sorted_result:
                if acc_num == self.__get_acc_prot_num(get_all_info(k[0])):
                    return k[0], k[1], k[3], True

        # if there is no length close enough to chain_length:
        if len(length_result) == 0:
            return "None", "None", "None", False

        # get the best identity from the length_result
        sorted_result = sorted(length_result, key=lambda x: x[1], reverse=True)
        i = 0
        names = []
        while i < len(sorted_result):
            if len(names) > 0 and identity != sorted_result[i][1]: break
            names.append(sorted_result[i][0])
            identity = sorted_result[i][1]
            i += 1

        # get the best resolution out of the same identities of proteins.
        # if there is no X-RAY structure, take the first NMR
        resolution = 100
        protein = names[0]
        for ch in range(0, i):
            result = describe_pdb(names[ch])
            if result['expMethod'] == 'X-RAY DIFFRACTION':
                try:
                    if float(resolution) > float(result['resolution']):
                        resolution = result['resolution']
                        protein = names[ch]
                except:
                    continue

        # get the identity and sequence of the chosen protein
        for elem in sorted_result:
            if elem[0] == protein:
                return protein, elem[1], elem[3], False

    '''
    =================
    Help Functions for for the programmer
    =================
    '''

    def __get_hit_params(self,hit):
        '''
        Find the seperate protein of a given chain index, with given conditions

        Parameters:
        ----------
        hit : string
            Represents a hit from the plast
        chain_length: int
            Represent the length of the chain of interest
        Returns
        -------
        pdb_name : string
            A 4 character string represent a PDB entry of the hit
        identity : int
            The identity percent of the hit to the chain of interest, calculated with BLAST
        length : int
            The length of the chain from the his
        seq : string
           The sequence of the hit of interest  
        '''
        for data in hit:
            if "pdbid" in data:
                if data.split(":")[1] != "1":
                    return "", "", "", ""
                pdb_name = data.split(":")[0]
                break
        for data in hit:
            if "Score" in data: break

        seq_list = data.split("Sbjct: ")
        seq = ""
        for i in range(1, len(seq_list)):
            seq += seq_list[i].split("\n")[0][4:seq_list[i].split("\n")[0].find(" ", 5)]
        data = str(data)
        start = data.find("Score")
        end = data.find("Positives")
        length = int(data[data.find("Length"):start].split()[2])
        identity = (data[start:end].split()[14][1:-3])
        return pdb_name, length, identity, seq

    def __get_acc_comp_num(self,i):
        '''Get the accession number of a chain in compolex
            Parameters:
            ----------
            i : int 
                The index of the present checked macro-molecule entity (chain) in the PDB entry
                
            Returns
            -------
                A string that represent the accession number that was required
            
            '''
        if "macroMolecule" in self.all_info["polymer"][i].keys():
            if type(self.all_info["polymer"][i]["macroMolecule"]) is list:
                return self.all_info["polymer"][i]["macroMolecule"][0]["accession"]["@id"]
            return self.all_info["polymer"][i]["macroMolecule"]["accession"]["@id"]
        return "-"

    def __get_acc_prot_num(self,all_info):
        '''Get the accession number of a chain in compolex
           Parameters:
           ----------
           all_info : string 
               The information about the protein,using "get_all_info()" function of pypdb

           Returns
           -------
               A string that represent the accession number that was required

           '''
        if "macroMolecule" in all_info["polymer"].keys():
            if type(all_info["polymer"]["macroMolecule"]) is list:
                return all_info["polymer"]["macroMolecule"][0]["accession"]["@id"]
            return all_info["polymer"]["macroMolecule"]["accession"]["@id"]
        return "-"


    def get_search_params(self,index):
        '''Get the parameters needed for the search of the seperate protein
        Parameters:
        ----------
        index : int 
            The index of the present checked chein in the get_all_info()["polymer"]
        Returns
        -------
        chain_id : string
            A string that contains the chains of a macro-molecule entity    
        chain_length : string
            A string that contains the chains of a macro-molecule entity    
        acc_num : string
            A string that contains the chains of a macro-molecule entity    
        '''
        chain_id=self.__get_chains(self.all_info["polymer"][index]["chain"])
        chain_length= int(self.all_info["polymer"][index]["@length"].strip())
        url_scop = "http://www.rcsb.org/pdb/explore/explore.do?structureId=" + self.pdb_id
        soup = BeautifulSoup(urlopen(url_scop).read(), "html.parser")
        mydivs = soup.find_all("div", class_="row")
        acc_num=mydivs[5+3*index].text.split("\n")[3].strip()
        #acc_num=self.__get_acc_comp_num(index)
        return chain_id,chain_length,acc_num

    def __get_chains(self,dict):
        '''Extracts the charachters represnt the chains of macro-molecule entity, given as dictionary
        Parameters:
        ----------
        dict : dictionary 
            A dictionary that contains the chains of a macro-molecule entity ad @id:char
        Returns
        -------
        final : string
            A string that contains the chains of a macro-molecule entity        
        '''

        final = ""
        if len(dict) == 1:
            final += dict["@id"]
        else:
            for i in dict:
                final += i["@id"]
        if type(final) is list:
            final = ''.join(final)
        return final

    def __sort(self,results, chain_length, acc_num):
        # sort the list according to length
        sorted_result = sorted(results, key=lambda x: x[2], reverse=True)

        # the identity is limited to 80, therefore the limit of length is 20% of the given length
        limit = int(chain_length / 10) * 2
        len_range = 0
        i = 0
        length_result = []
        flag = False
        # find the closest length to the chain_length
        while len(length_result) == 0 and len_range <= limit:
            i = 0
            while i < len(sorted_result) and sorted_result[i][2] > chain_length + len_range:
                i += 1
            # sorted_result[i][2]>=chain_length+len_range and
            while i < len(sorted_result) and sorted_result[i][2] <= chain_length + len_range and sorted_result[i][
                2] >= chain_length - len_range:
                length_result.append(sorted_result[i])
                i += 1
                if i == len(sorted_result) or not (sorted_result[i][2] <= chain_length + len_range and sorted_result[i][
                    2] >= chain_length - len_range):
                    flag = True

            if flag == True:
                break
            len_range += 1

        # no protein in the limited length
        if len(length_result) == 0:
            for k in sorted_result:
                if acc_num == get_acc_num(get_all_info(k[0])):
                    return [k], True

        # if there is no length close enough to chain_length:
        if len(length_result) == 0:
            return [], False

        # get the best identity from the length_result
        sorted_result = sorted(length_result, key=lambda x: x[1], reverse=True)
        return sorted_result, False

if __name__=="__main__":
    c=Complex("1c0t")
    print (c)
    pdbs=Search().get_list_from_file(r"C:\Users\Efrat\Desktop\mini project\Project\test.txt")
    shirit=Search().get_list_from_file(r"C:\Users\Efrat\Desktop\mini project\PDB\test.txt")
    print (len(pdbs))
    for i in shirit:
        if not i in pdbs:
            pdbs.append(i)

    print (len(pdbs),pdbs)

