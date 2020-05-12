import math
#mldist = "../data.tre.mldist"
oral_rep_file = "00.data/oral_sgb_representative.tsv"
tax_file = "00.data/oral_sgb_gtdb_taxonomy.tsv"
annot_file = "00.data/oral_annot.txt"
nodes_dict_file = "00.data/rename.dict"

nodes_dict = {}
for i in open(nodes_dict_file):
    k,v = i.strip().split("\t")
    nodes_dict[k] = v
print("nodes_rename_dict:", len(nodes_dict))

rep_dict = {}
for i in open(oral_rep_file).readlines()[1:]:
    item = i.split("\t")
    size, mtype, path = item[1], item[2], item[-2]
    rep_id = path.split("/")[-1].rstrip(".gz").rstrip(".fna").rstrip(".fa")
    rep_dict[rep_id] = [size, mtype]
print("load mgs_dict:", len(rep_dict))

tax_dict = {}
for i in open(tax_file).readlines()[1:]:
    item = i.split("\t")
    phylum = item[5].split(";")[1].split("p__")[1].split("_")[0]
    rep_id = item[9].split("/")[-1].rstrip(".gz").rstrip(".fna").rstrip(".fa")
    tax_dict[rep_id] = phylum
print("load tax_dict:", len(tax_dict))

f_o = open(annot_file, "w")
global_set = \
"total_plotted_degrees\t370\n\
start_rotation\t90\n\
ring_label\t1\tUnknownnes\n\
ring_label\t2\tGenomes\n\
ring_width\t2\t3\n\
ring_height\t1\t0.8\n\
ring_color\t2\t#949494\n\
ring_internal_separator_thickness\t1\t0.3\n\
ring_internal_separator_thickness\t2\t0.3\n\
clade_marker_size\t1.0\n\
branch_thickness\t0.5\n\
total_plotted_degrees\t340\n\
clade_marker_edge_color\t#555555\n\
clade_marker_edge_width\t0.1\n"
f_o.write(global_set)

#node_color
node_color = \
{"Firmicutes":"#CAB2D6",\
"Patescibacteria":"#FF7F00",\
"Actinobacteriota":"#FDBF6F",\
"Bacteroidota":"#E31A1C",\
"Proteobacteria":"#FB9A99",\

"Campylobacterota":"#33A02C",\
"Fusobacteriota":"#B2DF8A",\
"Spirochaetota":"#1F78B4",\
"Synergistota":"#A6CEE3",\
"Others":"#808080",\

"kSGB":"#4a87cb",\
"uSGB":"#e66827"}

# legend
for k,v in node_color.items():
    nstr_l = "%s\tannotation\t%s\n" %(k, k)
    nstr_l += "%s\tclade_marker_color\t%s\n" %(k, v)
    nstr_l += "%s\tclade_marker_size\t50\n" %(k)
    f_o.write(nstr_l)

# nodes
for node in nodes_dict.keys():
    phylum = tax_dict[node]
    if phylum in node_color.keys():
        pass
    else:
        phylum = "Others"
    nstr_node="%s\tclade_marker_color\t%s\n%s\tclade_marker_size\t10\n" %(nodes_dict[node], node_color[phylum], nodes_dict[node])
    nstr_node += "%s\tannotation_background_color\t%s\n" %(nodes_dict[node], node_color[phylum])
    nstr_node += "%s\tannotation_background_alpha\t0.1\n" %(nodes_dict[node])
    f_o.write(nstr_node)

#ring 1 unknownness
r1_color = {"kMGS":"#4a87cb", "uMGS":"#e66827"}
#r1_color = {"kSGB":"#AAAA00", "uSGB":"#AA00AA"}
for node in nodes_dict.keys():
    mtype = rep_dict[node][1]
    nstr_r1 = "%s\tring_color\t1\t%s\n" %(nodes_dict[node], r1_color[mtype])
    f_o.write(nstr_r1)

#ring 2 genomes(log10)
for node in nodes_dict.keys():
    size = math.log10(float(rep_dict[node][0]))/3
    nstr_r2 = "%s\tring_height\t2\t%.4f\n" %(nodes_dict[node], size)
    f_o.write(nstr_r2)
f_o.close()
