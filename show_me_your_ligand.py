# imports
import tempfile
import math
import subprocess
import os
import sys

import getpass
user = getpass.getuser()

try:
    from pymol import cmd
    import math
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)
try:
    cmd.set('cartoon_gap_cutoff', 0)
except:
    pass

print ("\n\n\n"+"Program create a list of residues interaction with hetero atoms (ligands) in radius of 5A:"+
       "\n"+"dirtyC - for quick clean the structures and superimpose them"
       +"\n"+"show_lig - marks ligands in molecule and make it visible :)"
       +"\n"+"show_AS - showing  residues in 5 A radius of ligand and coloring it on yellow. "
       +"\n"+"aa - superimpose all pdb files""\n\n")

def align_all( subset = [] ):
  """
  Superimpose all open models onto the first one.
  This may not work well with selections.
  This function is taken from rnapdbtool(:P) and before the function is probably taken form: https://daslab.stanford.edu/site_data/docs_pymol_rhiju.pdf
  """
  print("""This returns a list with 7 items:
    RMSD after refinement
    Number of aligned atoms after refinement
    Number of refinement cycles
    RMSD before refinement
    Number of aligned atoms before refinement
    Raw alignment score
    Number of residues aligned """)

  AllObj=cmd.get_names("all")
  for x in AllObj[1:]:
    #print(AllObj[0],x)
    subset_tag = ''
    if isinstance( subset, int ):
      subset_tag = ' and resi %d' % subset
    elif isinstance( subset, list ) and len( subset ) > 0:
      subset_tag = ' and resi %d' % (subset[0])
      for m in range( 1,len(subset)): subset_tag += '+%d' % subset[m]
    elif isinstance( subset, str ) and len( subset ) > 0:
      subset_tag = ' and %s' % subset
    values = cmd.align(x+subset_tag,AllObj[0]+subset_tag)
    print(AllObj[0], x, ' '.join([str(v) for v in values]), '-- RMSD', values[3], ' of ', values[6], 'residues')
    cmd.zoom()

def dirty_clean():
    """
    removing all chains (and water) beside chain A, adding hydrogens
    """
    print(' >>>>removing all chains (and water) beside chain A, adding hydrogens <<<<');
    cmd.do("remove all and not chain A");
    print(' >>>> removing all chains not A <<<<');
    cmd.do("select all");
    print(' >>>> adding Hydrogens <<<<');
    cmd.do("set_name sele, work_sele");
    cmd.do("h_add work_sele");
    cmd.do("h_add (work_sele)");
    cmd.do("remove resn hoh"); ## or remove solvent
    cmd.do("color gray")
    align_all(subset='all')


def show_ligands():
    """
    Select, shows ligands
    """

    print(' >>>> Show ligands <<<<');
    cmd.do('select HETATM');
    cmd.do('color blue, sele');
    cmd.do("set_name sele, ligands");
    cmd.do("util.cbac ligands");
    cmd.do("show spheres, ligands");
    cmd.do("set sphere_transparency, 0.6");
    cmd.do("zoom ligands, 6"); #this can be opptional , if you have few structures with ligands in many places

def show_active_site():
    """
    Select, shows active site of ligand in range - works for one pdb ...
    """
    print(' >>>> Show active site <<<<');
    cmd.do("AS_res_list = []")
    cmd.do("AS_L = []")
    #cmd.do("select active_site, ligands around 5");
    cmd.do("for n in cmd.get_names_of_type('object:molecule'): cmd.set_name(n,n[:5])")
    cmd.do("select active_site, (br. all within 5 of ligands) and not ligands")
    cmd.do("show sticks, active_site");
    cmd.do("util.cbay active_site");
    #cmd.do("iterate(active_site, ca), list.append((resi, resn))");
    cmd.do("iterate active_site, AS_res_list.append((resi, resn))");
    cmd.do("AS_L=(list(set(AS_res_list)))");
    cmd.do("for x in cmd.get_names_of_type('object:molecule'): cmd.set_name(x,x[:5])") # 5 A
    cmd.do("print (x,AS_L)");

def show_active_site1():
    """
    Select, shows active site of ligand in range of 5 Angstrom
    """
    print(' >>>> Show active site <<<<');
    file_object = open('Full_AS_intr_res_list.csv', 'w')
    file_object1 = open('AS_intr_res_list.csv', 'w')
    #cmd.do("select active_site, ligands around 5");
    obj_list = cmd.get_names('objects')
    print(obj_list)
    for obj in obj_list:
        cmd.do("AS_res_list=[]")
        cmd.do("AS_L = []")
        cmd.do("AS_L_num = []")
        cmd.do("Lig_list=[]")
        cmd.do("fix_Lig_list=[]")
        cmd.do("AS_res_num_list=[]")
        obj_us_list=[]
        print(obj)
        cmd.do("select HETATM and "+ str(obj))
        cmd.do("iterate sele, Lig_list.append((resi + resn))");
        #print list of ligands in pdbs
        cmd.do("fix_Lig_list=set(list(Lig_list))")
        #print((str(obj), fix_Lig_list))
        cmd.do("select sele1,( br. " + str(obj)+" within 5 of sele) and not sele")
        #and not sele")
        cmd.do("util.cbay sele1");
        cmd.do("iterate sele1, AS_res_list.append((resi, resn))");
        cmd.do("iterate sele1, AS_res_num_list.append((resi))");
        cmd.do("AS_L=(list(set(AS_res_list)))");
        cmd.do("AS_L_num=(list(set(AS_res_num_list)))");
        # print(str(obj), fix_Lig_list, fix_Lig_list, AS_L)\
        # change obj to list  for purpose of futher file playing :P
        cmd.do("if AS_L==[]: AS_L=['0']")
        cmd.do("if AS_L_num==[]: AS_L_num=['0']")
        file_object.write(csv_formatted_output(obj,fix_Lig_list, AS_L))
        file_object1.write(csv_formatted_output(obj,fix_Lig_list, AS_L_num))

        cmd.set_name("sele1", "AS_" + str(obj));

    file_object.close();
    file_object1.close();


def show_active_siteMYR():
    """
    Select, shows active site of MYR ligand in range of 5 Angstrom
    """
    print(' >>>> Show active site <<<<');
    file_object = open('Full_AS_intr_res_list.csv', 'w')
    file_object1 = open('AS_intr_res_list.csv', 'w')
    #cmd.do("select active_site, ligands around 5");
    obj_list = cmd.get_names('objects')
    print(obj_list)
    for obj in obj_list:
        cmd.do("AS_res_list=[]")
        cmd.do("AS_L = []")
        cmd.do("AS_L_num = []")
        cmd.do("Lig_list=[]")
        cmd.do("fix_Lig_list=[]")
        cmd.do("AS_res_num_list=[]")
        obj_us_list=[]
        print(obj)
        cmd.do("select resn MYR and "+ str(obj))
        cmd.do("iterate sele, Lig_list.append((resi + resn))");
        #print list of ligands in pdbs
        cmd.do("fix_Lig_list=set(list(Lig_list))")
        #print((str(obj), fix_Lig_list))
        cmd.do("select sele1,( br. " + str(obj)+" within 5 of sele) and not sele")
        #and not sele")
        cmd.do("util.cbay sele1");
        cmd.do("iterate sele1, AS_res_list.append((resi, resn))");
        cmd.do("iterate sele1, AS_res_num_list.append((resi))");
        cmd.do("AS_L=(list(set(AS_res_list)))");
        cmd.do("AS_L_num=(list(set(AS_res_num_list)))");
        #print(str(obj), fix_Lig_list, fix_Lig_list, AS_L)\
        #change obj to list  for purpose of futher file playing :P
        cmd.do("if AS_L==[]: AS_L=['0']")
        cmd.do("if AS_L_num==[]: AS_L_num=['0']")
        file_object.write(csv_formatted_output(obj,fix_Lig_list, AS_L))
        file_object1.write(csv_formatted_output(obj,fix_Lig_list, AS_L_num))

        cmd.set_name("sele1", "AS_" + str(obj));

    file_object.close();
    file_object1.close();


def csv_formatted_output(pdb_id, lig, actives):
    #makes pretty csv output file >pdbid >ligands >list of numbers: am. ac. of active site
    output = pdb_id+','
    output += ' '.join(lig) + ','
    if type(actives[0]) == tuple:
        actives = [i+' '+j for i, j in actives]
    output += ','.join(actives) + '\n'
    return output

def quickref():
    print('  dirtyC - fast clean of object, superposition')
    print(' show_lig - marking ligands in the pdb')
    print(' show_AS - coloring ac. site on yellow, generate file with list of am. in 5 A around ligands')

try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
else:
    cmd.extend('dirtyC', dirty_clean)
    cmd.extend('show_lig', show_ligands)
    cmd.extend('show_AS', show_active_site1)
    cmd.extend('show_myr', show_active_siteMYR)
    cmd.extend('aa', align_all)
