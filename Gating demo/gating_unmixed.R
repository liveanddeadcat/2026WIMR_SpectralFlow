setwd("D:/Spectral flow seminar and workshop/TRM example/Gating demo")
# load packages
library(tidyverse)
library(CytoExploreR)
library(flowCore)

# set the path and load files
fcs_dir  = "./fcs"

gs <- cyto_setup(fcs_dir,
                 gatingTemplate = "ols_gating.csv")

# build the data transformer
trans <- cyto_transformer_logicle(gs)

# transform the data
gs <- cyto_transform(gs, trans = trans)


# gating
cyto_gate_draw(gs,
               alias = "Time",
               channels = c("Time","SSC-A"),
               type = "rectangle",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "Time",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "Single Cells",
               alias = "AFlowCD45high",
               channels = c("AF","CD45"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "AFlowCD45high",
               alias = "Live",
               channels = c("LD","FSC-A"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "Live",
               alias = "CD3+HLADRlow",
               channels = c("CD3","HLADR"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")


cyto_gate_draw(gs,
               parent = "CD3+HLADRlow",
               alias = "T cells",
               channels = c("FSC-A","SSC-A"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")


cyto_gate_draw(gs,
               parent = "T cells",
               alias = "CD4+ T cells",
               channels = c("CD4","CD8"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "CD4+ T cells",
               alias = "CD69+",
               channels = c("CD69", "CD103"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")

cyto_gate_draw(gs,
               parent = "CD4+ T cells",
               alias = "CD103+",
               channels = c("CD69", "CD103"),
               type = "polygon",
               gatingTemplate = "ols_gating.csv")


# gating strategy
cyto_plot_gating_scheme(gs, label_text_size = 0.6)

# Scatter plot
cyto_plot(gs, 
          parent = "CD4+ T cells",
          alias = "CD69+",
          channels = c("CD69","CD103"),
          border_line_width = 0.7,
          axes_text_font = 1,
          xlim = c(-1e4, 1e6),
          ylim= c(-1e4, 1e6))


# Directly extract gated CD4+ T cell logicle-transformed values
cd4_t <- cyto_extract(gs, "CD4+ T cells", raw = TRUE)

# Save raw values
cd4_t <- cyto_extract(gs, "CD4+ T cells")
cd4_t_ff <- cyto_convert(gs, return = "flowFrame", parent ="CD4+ T cells") # convert to flowframe
cd4_t_ff_raw <- cyto_transform(cd4_t_ff, trans = trans,inverse = TRUE) # inverse tranform
  
#expression data
cd4_df_raw <- exprs(cd4_t_ff_raw)
summary(cd4_df_raw)

# Save the gating set object
cyto_save(gs, save_as = "ols_gs_colonLP")
# read the saved gs.
gs <- cyto_load(path="ols_gs_colonLP")

# Summary statistics
cyto_stats_compute(gs,
                   alias = c("CD69+","CD103+"),
                   parent = "CD4+ T cells",
                   stat = "freq")
cyto_stats_compute(gs,
                  alias = c("CD69+","CD103+"),
                  channels = c("CD69"),
                  stat = "geo mean")
cyto_stats_compute(gs,
                   alias = c("CD69+","CD103+"),
                   channels = c("CD69"),
                   stat = "median")

