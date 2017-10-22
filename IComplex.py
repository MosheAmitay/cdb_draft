from pypdb.pypdb import *

class IComplex:

    """resolution =0
    pub_year =""
    pdb_id=""
    all_info=""
    description="""""



    def __init__(self,name_pdb):
        '''Initiating a Complex object, with given PDB code
           Parameters
           ----------
           name_pdb : string
               A 4 character string giving a pdb entry of interest

           '''
        self.pdb_id=name_pdb
        self.all_info=get_all_info(self.pdb_id)
        self.description=describe_pdb(self.pdb_id)
        self.resolution=self.set_resolution()
        self.pub_year=self.description['release_date'][0:4]

    def set_resolution(self,description=""):
        '''Set the resolution of complex, if possible

          Parameters:
          ----------
          description : string
              Get a description of proteien.        
              Get description="" as default, to use self.description
          Returns
          -------
          resolution : string
              A string represents the resolution of given PDB entry

         '''
        resolution=0
        if description=="":
            description=self.description
        try:
            if description['expMethod'] == 'X-RAY DIFFRACTION':
                resolution = description['resolution']
            return resolution
        except:
            return resolution


    def __str__(self):
        '''Print function
          Examples
          --------
          c=IComplex("1A0O")
          print (c)
          >>>The Complex: 1A0O (Publication Year: 1998, Resolution: 2.95)
        '''
        details = "The Complex: " + self.pdb_id + " (" + "Publication Year: " + self.pub_year + ", Resolution: " + self.resolution + ")"
        return details