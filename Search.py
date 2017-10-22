from pypdb.pypdb import *
from urllib.request import urlopen
from bs4 import BeautifulSoup

class Search():
    '''
   =================
   Help Functions for for the programmer
   =================
   '''

    def __findScop(self,pdb_id):
        '''
        Find if SCOP anotation of PDB entry contains Domain
        Parameters:
        ----------
        pdb_id : string 
              A 4 character string giving a pdb entry of interest

        Returns
        -------
            bool
            A boolean codition,indictes if there is "Domain" or not
       '''
        family = ""
        url_scop = "http://www.rcsb.org/pdb/explore/macroMoleculeData.do?structureId=" + pdb_id
        soup = BeautifulSoup(urlopen(url_scop).read(), "html.parser")
        mydivs = soup.find_all("div", class_="row")
        div = mydivs[5]
        if len(re.findall("Domain", mydivs[6].div.text)) >= 1:
            return True
        return False


    def __findScop_cath(self, pdb_id):
        '''
       Find if SCOP anotation of PDB entry contains Domain
       Parameters:
       ----------
       pdb_id : string 
             A 4 character string giving a pdb entry of interest

       Returns
       -------
           bool
           A boolean expression,indictes if there is "Immunoglobulin-like" in cath annotation, or not
      '''
        family = ""
        url_scop = "http://www.rcsb.org/pdb/explore/macroMoleculeData.do?structureId=" + pdb_id
        soup = BeautifulSoup(urlopen(url_scop).read(), "html.parser")
        mydivs = soup.find_all("div", class_="row")
        div = mydivs[5]
        if div.text.find("Immunoglobulin-like") != -1:
            return True
        return False

    def __checkWords(self,pdb_id,direcotry):
        '''
        Check header words
        Parameters:
       ----------
       pdb_id : string 
             A 4 character string giving a pdb entry of interest
        directory : string
            Where to keep the downloaded structure

       Returns
       -------
        bool
           A boolean expression,indictes if there is "antibody","fab" or "immunoglob" or not
        '''
        url_scop = "https://files.rcsb.org/view/"+pdb_id+".pdb"
        soup = BeautifulSoup(urlopen(url_scop).read(), "html.parser")
        soup = soup.text.split("KEYWDS")
        text=""
        for i in range(1, len(soup)):
            text+=soup[i].split("\n")[0]+" "
        text=text.lower()
        if text.find("antibody") != -1 or text.find("fab") != -1 or text.find("immunoglob") != -1:
            return True
        return False

        '''
     =================
     Functions for for the user
     =================
     '''
    def search_complexes(self, pdbs, num_chains=2, length_chains=30, title_flag=True, counter=0, counter_end=0):
        '''Get the complexes from a list of proteins (pdb codes), with given conditions and range

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
           '''

        complexes = []
        errors = []

        if counter_end == 0:
            counter_end = len(pdbs)
        for i in range(counter, counter_end):
            pdb = pdbs[i][0:4]
            all_info=get_all_info(pdb)
            description=describe_pdb(pdb)
            counter += 1
            try:
                if (len(all_info["polymer"]) == num_chains):
                    flag = True
                    for i in range(num_chains):
                        this_length = int(all_info["polymer"][i]["@length"].strip())
                        if this_length < length_chains:
                            flag = False
                            break
                    if flag == False: continue
                    all_names = []
                    for i in range(num_chains):
                        all_names.append(all_info["polymer"][i]["polymerDescription"]["@description"])

                    for t in range(len(all_names) - 1):
                        for j in range(t+1, len(all_names)):
                            if all_names[t] == all_names[j]:
                                t = len(all_names)
                                flag = False
                                break

                    if flag == False: continue
                    title = description["title"].lower()
                    if title_flag == True:
                        if not ("complex" in title) and not ("bound" in title):
                            continue
                    else:
                        if prot1_name.upper().find("PROT") == -1 or prot2_name.upper().find("PROT") == -1:
                            try:
                                if self.__findScop(pdb) == False:
                                    continue
                            except:
                                print("Error occures in finding the SCOP. Try check this protein again. Protein code: " + self.pdb_id)

                    complexes.append(pdb)
            except:
                errors.append(pdb)
        return complexes, errors

    def get_list_from_file(self, file_name=""):
        '''Extracts a list from a file of pdb entries
         Parameters:
         ----------
         file_name : string 
             A string that represnt file name
             Gets file_name="" as default in order to be replaced with get_all() in the function,
             which is all the proteins in the PDB
         Returns
         -------
         pdbs : list
             A list where each value is a 4 character string giving a pdb entry of interest
        '''
        if file_name == '':
            return get_all()
        file = open(file_name, 'r')
        pdbs = file.read().split('\n')
        new=[]
        file.close()
        #i = len(pdbs)-1
        #delete empty cells
        for i in pdbs:
            if len(i)>0:
                new.append(i)
        return pdbs


    def search_antibodies(self, pdbs, num_chains_min=2, num_chains_max=6,length_chains=30, counter=0, counter_end=0):
        '''Get the antibodies from a list of proteins (pdb codes), with given conditions and range

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
        '''
        antibodies = []
        errors = []
        if counter_end == 0:
            counter_end = len(pdbs)
        for i in range(counter, counter_end):
            pdb = pdbs[i]
            all_info = get_all_info(pdbs[i])
            counter += 1
            #try:
            """check how many entities in macromolcule. if there ara more than two:"""
            try:
                if not type(all_info["polymer"]) is list:
                    continue

                if len(all_info["polymer"]) > num_chains_min and len(all_info["polymer"]) <= num_chains_max:
                    description = []

                    """check if classification is antibody"""
                    if describe_pdb(pdbs[i])['keywords'].upper().find("ANTIBODY") != -1:
                        antibodies.append(pdb)


                    else:
                        k = 0
                        for j in range(len(all_info["polymer"])):
                            if int(all_info["polymer"][j]["@length"].strip())<length_chains:
                                k = len(all_info["polymer"])
                                break
                        """check if polymer descriptions contains antibody"""
                        for j in range(len(all_info["polymer"])):
                            this_des = all_info["polymer"][j]["polymerDescription"]["@description"]
                            description.append(all_info["polymer"][j]["polymerDescription"]["@description"])
                            if this_des.upper().find("ANTIBODY") != -1:
                                antibodies.append(pdb)
                                k = len(all_info["polymer"])
                                break

                        """if not-check if macromolecule description are similar, and if so-check cath classification"""
                        y = k + 1
                        while k < len(all_info["polymer"]) - 1:
                            while y < len(all_info["polymer"]):
                                if description[k].find(description[y][0:3]) != -1 or description[k].find(
                                        description[y][len(description[y]) - 3:len(description[y])]) != -1:
                                    if self.__findScop_cath(pdbs[i]) == True:
                                        directory=r"C:\Users\Efrat\Desktop\mini project\Project"
                                        if self.__checkWords(pdbs[i],directory):
                                            antibodies.append(pdb)
                                        k = y = all_info["polymer"]
                                        break
                                    else:
                                        k = y = all_info["polymer"]
                                        break
                                y += 1
                            if k == y: break
                            k += 1
                            y = k + 1
            except:  # if No Annotations Available
                errors.append(pdb)
        return antibodies,errors






