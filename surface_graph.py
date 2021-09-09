
class Surface_Amino_Acid:

    def __init__(self, AA_name, index):
        self.AA_name = AA_name
        self.index = index
        self.Neighbors_list = []
        self.part_of_the_path = False


    def add_new_neighbor(self , neighbor):
        self.Neighbors_list.append(neighbor)


    def insert_to_path(self):
        self.part_of_the_path = True

    def check_if_in_path(self):
        if self.part_of_the_path == True:
            return (True)

    def get_neighbors_list(self):
        lst=[]
        for neighbor in self.Neighbors_list:
            lst.append(neighbor.AA_name + neighbor.index)
        return(lst)
            #return(neighbor.AA_name + neighbor.index)
