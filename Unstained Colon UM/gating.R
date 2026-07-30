# load packages
library(tidyverse)
library(CytoExploreR)
library(flowCore)

# set the path and load files
fcs_dir  = "./fcs"

gs <- cyto_setup(fcs_dir,
                 gatingTemplate = "unstained_gate.csv")

# build the data transformer
trans <- cyto_transformer_logicle(gs)

# transform the data
gs <- cyto_transform(gs, trans = trans)


# gating
cyto_gate_draw(gs,
               alias = "Time",
               channels = c("Time","SSC-A"),
               type = "rectangle",
               gatingTemplate = "unstained_gate.csv")

cyto_gate_draw(gs,
               alias = "Cells",
               parent = "Time",
               channels = c("FSC-A","SSC-A"),
               type = "polygon",
               gatingTemplate = "unstained_gate.csv")


## export files

#AF profile
af <- cyto_extract(gs, "Cells")
af_ff <- cyto_convert(af, return = "flowFrame",
                          inverse.transform = TRUE)


# Output fcs files
write.FCS(x = af_ff, filename = "Lympho_af.fcs")



