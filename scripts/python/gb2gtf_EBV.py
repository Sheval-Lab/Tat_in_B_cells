#!/usr/bin/env python
# coding: utf-8

# This notebook is used to convert EBV gene annotation file in the GenBank format into GTF format. It does not retain all annotated features, it keeps only annotated coordinates of exons and genes (if exons coordinates are missing). The notebook produces GTF file that can be used for EBV gene expression counting (with `type = exon`).

# In[1]:


from Bio import SeqIO


# In[2]:


# GenBank file as input annotation file
input_gb = "../../data/metadata/EBV.gb"

# GTF file to write output
output_gtf = "../../data/metadata/EBV_exons.gtf"


# In[3]:


with open(input_gb) as handle:
    for record in SeqIO.parse(handle, "gb"):
        # Extract 'chromosome' name
        seqname = record.id
        
        # allow to process only required types: genes and exons
        types2process = ["gene", "mRNA"]
        
        # Create an empty list to save each processed feature
        gtf = []
        
        # Create an empty list to keep track of genes with annotated exons
        exons = []
        
        for f in record.features:
            # Process only specified feature types
            if f.type not in types2process:
                continue
            
            # Skip feature if it does not contain gene_name identifier
            if "gene" not in f.qualifiers:
                continue
                
                
            # If feature is OK, process it
            # Extract gene_name
            gene_name = f.qualifiers['gene'][0]
            # Extract strand as + or -
            strand = "+" if f.strand == 1 else "-"
            # Extract feature type: gene or exon
            f_type = f.type
            # Extract feature coordinates
            coords = f.location.parts
            
            # Save gene_name if this gene has annotated exons
            if f_type == "mRNA":
                exons.append(gene_name)
            
            # For each exon or gene create an individual line in GTF file
            for part in coords:
                gtf.append({"seqname": seqname, "gene_name": gene_name ,"type": f_type, "start": int(part.start), "end": int(part.end), "strand": strand, "index": coords.index(part)})                
  


# In[4]:


source = "gb2gtf"
f_type = "exon"

# Write GTF lines into file: keep only exons and genes without annotated exons
with open(output_gtf, 'a') as out_gtf:
    for record in gtf:
        if record["gene_name"] not in exons:
            gtf2file = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % (record["seqname"], source, f_type, record["start"]+1, record["end"], record["strand"], 'gene_id "'+record["gene_name"]+'"; exon_id "'+record["gene_name"]+'_'+str(record["index"])+'";')
            out_gtf.write(gtf2file)
        if record["gene_name"] in exons and record["type"] == "mRNA":
            gtf2file = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % (record["seqname"], source, f_type, record["start"]+1, record["end"], record["strand"], 'gene_id "'+record["gene_name"]+'"; exon_id "'+record["gene_name"]+'_'+str(record["index"])+'";')
            out_gtf.write(gtf2file)
        else:
            continue
        

