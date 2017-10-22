from pypdb.pypdb import *
from IComplex import IComplex
from Search import *
from Complex import Complex
from urllib.request import urlopen
from bs4 import BeautifulSoup


class Antibody(Complex):
    '''
    =================
    Functions for for the user
    =================
    '''

    def __init__(self, name_pdb):
        '''Initiating a Complex object, with given PDB code. Using the Complex initiation.
           Parameters
           ----------
           name_pdb : string
               A 4 character string giving a pdb entry of interest

           '''
        Complex.__init__(self, name_pdb)
        heavy,light,proteins=self.get_indexes()
        self.__heavy=heavy
        self.__light=light
        self.__proteins=proteins
        self.__sep_prot_details={}

    def get_indexes(self):
        '''Get the indexes of the antibodies, and the seperate protein
           Returns
           -------
           light : list
                A list with the indexes of the light chains in the PDB entry, in respective to the heavy chain
           heavy : list
                A list with the indexes of the heavy chains in the PDB entry, in respective to the light chain
           proteins : list
                A list with the indexes of the chains of the seperate protein of the PDB entry
           Examples
           --------
            a=Antibody("1P2C")
            print (a.get_indexes())
           >>>([0], [1], [2])
           '''
        description = []
        for j in range(len(self.all_info["polymer"])):
            description.append(self.all_info["polymer"][j]["polymerDescription"]["@description"])
        if len(self.all_info["polymer"]) > 4:
            light, heavy, proteins=self.__find_index_for_many_antibody(description)
        else:
            light, heavy, proteins= self.__find_index_for_one_antibody(description)
        return light,heavy,proteins



    def find_seperate_proteins_anti(self):
        '''Find the seperate proteins of antibody

        Returns
        -------
        prot_names : list
            A list of PDB 4 charachters entry, each represents seperate protein of chain.

        Examples
        --------
        b=Antibody("3TVM")
        b.find_seperate_proteins_anti()
        >>>['4ZAK', '3TVM']
        '''
        prot_names = []
        # for the seperate protein
        if len(self.__proteins)==1:
            pdb_code, identity, seq, flag = self.find_protein(self.__proteins[0])
        else:
            pdb_code, identity, seq, flag = self.find_protein_anti(self.__proteins[0],self.__proteins[1])
        if pdb_code == "None":
            self.__sep_prot_details = {}
            self.__sep_prot_details["None"] = "none"
            return []
        this_details = (identity, seq, flag)
        self.__sep_prot_details[pdb_code] = this_details
        prot_names.append(pdb_code)
        for i in range(len(self.__heavy)):
            #for the complex
            pdb_code, identity, seq, flag = self.find_protein_anti(self.__heavy[i],self.__light[i])
            if pdb_code == "None":
                self.__sep_prot_details = {}
                self.__sep_prot_details["None"] = "none"
                return []
            this_details = (identity, seq, flag)
            self.__sep_prot_details[pdb_code] = this_details
            prot_names.append(pdb_code)
        return prot_names

    def find_protein_anti(self,chain1,chain2,identity_percent=80,limit_percent=20):
        '''
        Find the seperate protein of the antibody's complexed protein, in case this protein itself is complex
        Parameters:
        ----------
        chain1: string
            The first chain of the complexed protein
        chain2: string
            The second chain of the complexed protein       
        identity_percent: int
             The identity percent that the result will be filteres accordin to. 
             Gets 80% as Default
        limit_percent: int
            The limit percent for cahin length that the result will be filteres accordin to.
             Gets 20% as different

        Returns
        -------
            The result of calling the function find_protein_many,which is name, identity, sequence and flag of the complexed protein's separate protein
        '''
        chain_id1, chain_length1, acc_num1 = self.get_search_params(chain1)
        chain_id2, chain_length2, acc_num2 = self.get_search_params(chain2)
        return self.find_protein_many(chain_id1, chain_id2, chain_length1, chain_length2, acc_num1, acc_num2,identity_percent,limit_percent)


    def find_protein_many(self, chain_id1, chain_id2, chain1_length1, chain1_length2, ac1, ac2, identity_percent, limit_percent):
        '''Find seperate protein of two chains. Example: light+heavy chain of antibody
        Paraameters:
        -----------
        chain_id1: string
            The first chain of the set of chains
        chain_id2: string
            The second chain of the set of chains
        chain1_length1: int
            The length of the first chain
        chain1_length2: int
            The length of the second chain
        ac1: string
            The accession number of the first chain
        ac2: string
            The accession number of the second chain
        identity_percent: int
             The identity percent that the result will be filteres accordin to. 
        limit_percent: int
            The limit percent for cahin length that the result will be filteres accordin to. 
        Returns:
        --------
        string
            4 charachters string, represents PDB entry of the seperate protein
        int     
            The identity of the separate proteins vs the checked 2 chains
        string
            The sequence of the separate protein
        bool
            To indicates wether the separate protein was found by accession number or not 
        Example:
        -------
        b=Antibody("3TVM")
        b.find_protein_many('A','B',285,99, "P11609","P01887",80,20)
        >>>('4ZAK', 100.0, 'SEAQQKNYTFRCLQMSSFANRSWSRTDSVVWLGDLQTHRWSNDSATISFTKPWSQGKLSNQQWEKLQHMFQVYRVSFTRDIQELVKMMSPKEDYPIEIQLSAGCEMYPGNASESFLHVAFQGKYVVRFWGTSWQTVPGAPSWLDLPIKVLNADQGTSATVQMLLNDTCPLFVRGLLEAGKSDLEKQEKPVAWLSSVPSSAHGHRQLVCHVSGFYPKPVWVMWMRGDQEQQGTHRGDFLPNADETWYLQATLDVEAGEEAGLACRVKHSSLGGQDIILYWQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWFYILAHTEFTPTETDTYACRVKHASMAEPKTVYWDRDM', False)
        '''
        # get blast of the 2 chains, separately
        result_blast1, self1 = self.__get_blast_result(chain_id1, chain_id2, identity_percent)
        result_blast2, self2 = self.__get_blast_result(chain_id2, chain_id1, identity_percent)

        # if there is no blast result
        if len(result_blast1) == 0 or len(result_blast2) == 0:
            return "None", "None", "None", False

        # if there is only self complex-return it
        if len(result_blast2) == 1 or len(result_blast1) == 1:
            if len(result_blast1) == 1:
                for j in range(len(result_blast2)):
                    if result_blast2[j][0] == result_blast1[0][0]:
                        return result_blast1[0][0], float(result_blast1[0][1] + result_blast2[j][1]) / 2, \
                               result_blast1[0][3] + result_blast2[j][3], False

                # return itself if not found
                if self1 == "" and self2[0] == pdb_id:
                    return self2[0], self2[1], self2[3], False
                elif self2 == "" and self1[0] == pdb_id:
                    return self1[0], self1[1], self1[3], False
                elif self1[0] == self.pdb_id and self2[0] == self.pdb_id:
                    return self1[0], self1[1], self2[3] + self1[3], False
                else:
                    return "None", "None", "None", False

            if len(result_blast2) == 1:
                for j in range(len(result_blast1)):
                    if result_blast1[j][0] == result_blast2[0][0]:
                        return result_blast2[0][0], float(result_blast2[0][1] + result_blast1[j][1]) / 2, \
                               result_blast2[0][3] + result_blast1[j][3], False
                # יש לך פה טעות נראלי!!! זה קוד שחזור על עצמו!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if self1 == "" and self2[0] == pdb_id:
                    return self2[0], self2[1], self2[3], False
                elif self2 == "" and self1[0] == pdb_id:
                    return self1[0], self1[1], self1[3], False
                elif self1[0] == self.pdb_id and self2[0] == self.pdb_id:  #############################הוספת self.
                    return self1[0], self1[1], self2[3] + self1[3], False
                else:
                    return "None", "None", "None", False
        # get all names that in 2 of the lists
        result1 = []
        result2 = []
        for name1 in result_blast1:
            for name2 in result_blast2:
                if name1[0] == name2[0]:
                    result1.append(name1)
                    result2.append(name2)
                    break

        if len(result1) == 0:
            if self1 == "" and self2 == "":
                return "None", "None", "None", False
            if self1 == "":
                if self2[0] == self.pdb_id:
                    return self2[0], self2[1], self2[3], False
            elif self2 == "":
                if self1[0] == self.pdb_id:
                    return self1[0], self1[1], self1[3], False
            elif self1[0] == self.pdb_id and self2[0] == self.pdb_id:
                return self1[0], self1[1], self1[3] + self2[3], False
            else:
                return "None", "None", "None", False

        sorted1, flag1 = self.__sort(result1, chain1_length1, ac1, limit_percent)
        sorted2, flag2 = self.__sort(result2, chain1_length2, ac2, limit_percent)

        for n1 in sorted1:
            for n2 in sorted2:
                if n1[0] == n2[0]:
                    if flag1 == True or flag2 == True:
                        return n1[0], float(n1[1] + n2[1]) / 2, n1[3] + n2[3], True
                    return n1[0], float(n1[1] + n2[1]) / 2, n1[3] + n2[3], False

        return result1[0][0], float(result1[0][1] + result2[0][1]) / 2, result1[0][3] + result2[0][3], False

    def get_proteins_details(self):
        '''Returns the details of the  seperate proteins of complex, refers to the complex

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
           '''
        if len(self.__sep_prot_details.keys()) == 0:  # the separated proteins wasn't found yet
            not_to_use = self.find_seperate_proteins()
        return self.__sep_prot_details

    def __str__(self):
        '''Print function
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

           '''
        if len(self.__sep_prot_details.keys()) == 0:  # the separated proteins wasn't found yet
            not_to_use = self.find_seperate_proteins_anti()
        sep_details=self.__sep_prot_details
        if "None" in self.__sep_prot_details.keys():
            return "There is no seperate protein for this complex"
        details_print = "%%%%%%" + IComplex.__str__(self) + "%%%%%%\n"
        protein_details=self.__get_protein_details_print()
        index = 0
        keys=list(sep_details.keys())
        for i in range(len(self.__heavy)):
            details_print+=self.__help_print("Protein",keys[0],sep_details[keys[0]],protein_details)
            index=[self.__light[i], self.__heavy[i]]
            details=self.get_antibody_details(index,keys[i+1])
            details_print+=self.__help_print("Antibody",keys[i+1],sep_details[keys[i+1]],details)

        return details_print

    def get_antibody_details(self, index, protein):
        '''Get the complexes deatails

        Parameters:
        ----------
        index : list
          The indexes of the chain of interest
        protein: string
          A 4 charachters string, represents the PDB entry for the seperate protein
        Returns
        -------
        details : list
          A list wit protein's details: name, accession numer,chains were checked, chains length,resolution,publication year
        '''
        details = []
        description = describe_pdb(protein)
        chain_all = ""
        name_all = ""
        acc_num_all = ""
        chain_length_all = 0
        for i in index:
            chain_id, chain_length, acc_num = self.get_search_params(i)
            chain_all += chain_id
            name_all = name_all + "," + self.all_info["polymer"][i]["polymerDescription"]["@description"]
            chain_length_all += chain_length
            acc_num_all = acc_num_all + "," + acc_num
        c = Complex(protein)
        resolution = c.resolution
        pub_year = description['release_date'][0:4]
        details.append(name_all)
        details.append(acc_num_all)
        details.append(chain_all)
        details.append(chain_length_all)
        details.append(resolution)
        details.append(pub_year)
        return details

    '''
      =================
      Help Functions for for the programmer
      =================
      '''

    def __similar(self, description, all_heavy, index):
        '''Get the similariest description, as light and heavy couple

           Parameters:
           ----------
           description : list
                All of the descriptions of macro-molecule entity
            all_heavy : list
                All of the heavy chains indexes in the antibody's PDB entry
           index: int
               The index of the heavy chain of interest
           Returns
           -------
           light : int
                The light index, that respectively matches the heavy index
           '''
        k = 0
        light = 0
        my_similar = -1
        while k < len(self.all_info["polymer"]):
            if not k in all_heavy:
                w1 = description[k].upper()
                w2 = description[all_heavy[index]].upper()
                w1 = w1 + ' ' * (len(w2) - len(w1))
                w2 = w2 + ' ' * (len(w1) - len(w2))
                similar = sum(1 if i == j else 0 for i, j in zip(w1, w2)) / float(len(w1))
                if similar > my_similar:
                    my_similar = similar
                    light = k
            k = k + 1
        return light

    def __find_index_for_one_antibody(self, description):
        '''Find the indexes of light, heavy lists, in case there is one pair of light-heavy chain in PDB entry

           Parameters:
           ----------
           description : list
                All of the descriptions of macro-molecule entity
           Returns
           -------
           light : list
                The light index, that respectively matches the heavy index
            heavy : list
                The heavy index, that respectively matches the light index
            proteins : list
                The index of the seperate protein
           '''
        # find light,heavy chain according to similarity in name of macro-molecule entity
        k = 0

        y = k + 1
        flag = False
        while k < len(self.all_info["polymer"]) - 1:
            while y < len(self.all_info["polymer"]):
                if description[k].find(description[y][0:3]) != -1 or description[k].find(
                        description[y][len(description[y]) - 3:len(description[y])]) != -1:
                    flag = True
                    break

                y += 1
            if flag: break
            k += 1
            y = k + 1
        # light and heavy chains were not found by similarity
        if flag == False:
            y = 0
            while y < len(all_info["polymer"]):
                if description[y].upper().find("LIGHT") != -1:
                    break
                y = y + 1
            k = 0
            while k < len(all_info["polymer"]):
                if description[k].upper().find("HEAVY") != -1:
                    break
                k = k + 1

        proteins = []
        for i in range(len(self.all_info["polymer"])):
            if i != k and i != y:
                proteins.append(i)
        heavy = [y]
        light = [k]
        return light, heavy, proteins

    def __find_index_for_many_antibody(self, description):
        '''Find the indexes of light, heavy lists, in case there is one pair of light-heavy chain in PDB entry

         Parameters:
         ----------
         description : list
              All of the descriptions of macro-molecule entity
         Returns
         -------
         light : list
              The light indexes, that respectively matches the heavy indexes
          heavy : list
              The heavy indexes, that respectively matches the light indexes
          proteins : list
              The index/es of the seperate protein
         '''

        num_anti = 0
        k = 0
        heavy = []
        while k < len(self.all_info["polymer"]):
            if description[k].upper().find("HEAVY") != -1:
                num_anti = num_anti + 1
                heavy.append(k)
            k = k + 1

        light = []
        for i in range(num_anti):
            light.append(self.__similar(description, heavy, i))
        proteins = []
        for i in range(len(self.all_info["polymer"])):
            if not i in heavy and not i in light:
                proteins.append(i)
        return light, heavy, proteins

    def __sort(self, results, chain_length, acc_num, limit_percent):
        '''
        Sort the blast result
        Parameters:
        -----------
        results: list
            List of tuples, with blast result, each tuple contains: name, identity, sequence, flag
        chain_length: int
            Indicates the length of the chain of interest
        acc_num: string
            The accession number of the chain of interest
        limit_percent: int
            The maximal  percent, to the range of chain length to filter with
        Returns:
        --------
        sorted_result: list
            The results in sorted way
         bool
            Boolean to indicates whether the results were found according to accession number 

        '''
        # sort the list according to length
        sorted_result = sorted(results, key=lambda x: x[2], reverse=True)

        # the identity is limited to 80, therefore the limit of length is 20% of the given length
        limit = int(chain_length / 10) * limit_percent
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

    def __get_hit_params_anti(self, hit):
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
                if int(data.split(":")[1]) > 2 and data.split(":")[0] != self.pdb_id:
                    return -1, -1, -1, -1
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

    def __get_blast_result(self, chain_id, chain2_id, identity_percent):
        '''
        Get the blast result of a chain, and initiate the myself object if possible
        Parameters:
        -----------
        chain_id: string
            The chain of interest
        chain2_id: string
            The other chain of the complex.Used to initiate the myself object
        identity_percent: int
            The ideneity percent to filter the result
        Returns:
        --------
        results: list
            List of tuples, each indicats blast result,and contains: name of blast result, identity to chain, length and sequence
        myself: tuple
            A tuple of the self chain, contains: name of blast result, identity to chain, length and sequence
        '''
        myself = ''
        # do blast
        blast_results = get_blast2(self.pdb_id, chain_id=chain_id)
        if len(blast_results) == 0:
            return []
        i = 0
        results = []
        # get all pdb_id with more than 80% identity
        for res in range(len(blast_results[1])):
            this_result = self.__get_hit_params_anti(blast_results[1][res])
            if int(this_result[2]) == -1:
                continue
            if this_result[0] == self.pdb_id:
                myself = (this_result[0], int(this_result[2]), int(this_result[1]), this_result[3])
                continue
            if int(this_result[2]) >= identity_percent:
                # name,identity,length,seq
                tup1 = (this_result[0], int(this_result[2]), int(this_result[1]), this_result[3])
                results.append(tup1)

        # no pdb_id with more than identity_percent,take yourself as seperate protein
        tup = ()
        if myself != '':
            if myself[0] == self.pdb_id:
                url_scop = "http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=" + myself[0]
                seq = ""
                try:
                    soup = BeautifulSoup(urlopen(url_scop).read(), "html.parser")
                    for text in soup.text.split(">"):
                        if len(text) == 0: continue
                        if text[5] == chain_id or text[5] == chain2_id:
                            seq += text[text.find("\n"):len(text)].replace('\n', '')

                except:
                    pass
                tup = (myself[0], myself[1], myself[2], seq)
                if len(results) == 0:
                    results.append(tup)

        return results, myself

    def __get_protein_details_print(self):
        '''
        Get protein details for print propose
        Returns
        -------
        String with the details for the printing
        '''
        if len(self.__proteins)==1:
            c=Complex(self.pdb_id)
            return c.get_complex_details(self.__proteins[0],list(self.__sep_prot_details.keys())[0])
        else:
            index=[]
            for i in self.__proteins:
                index.append(i)
            mykey=""
            for key in self.__sep_prot_details:
                mykey=key
                break
            return self.get_antibody_details(index, mykey)


    def __help_print(self,entity,key,value,details):
        '''
        Help the __str__function
        Parameters:
        ----------
        entity: string
            to print "Protein" or "Antibody" for each line    
        key: string
            A 4 charachters string, represents the 
        value: tuple
            Contains the key details (sequence,identity and flag)
        details:
            All the details about the proteins from get_antibody_details or __get_protein_details_print 
        Returns:
        -------
        details_print: string
            Contains the string for the __str__ function
        '''
        details_print=""
        newkey = key
        if value[2] == True: newkey += "*"
        details_print += "\n%%% "+entity+": " + newkey + " %%%"
        details_print += "\n-"+entity+" Name: " + details[0]
        details_print += "\n-"+entity+" Accession Number: " + details[1]
        details_print += "\n-"+entity+" Chains: " + details[2]
        details_print += "\n-"+entity+" Chain Length: " + str(details[3])
        details_print += "\n-"+entity+" Sequence: " + value[1]
        details_print += "\n-"+entity+" Identity to Complex's Chain/s: " + str(value[0])
        details_print += "\n-"+entity+" Resolution: " + str(details[4])
        details_print += "\n-"+entity+" Publication Year: " + details[5]+"\n"
        return details_print

if __name__=="__main__":
    a=Antibody("3TVM")
    print (a.get_antibody_details())
