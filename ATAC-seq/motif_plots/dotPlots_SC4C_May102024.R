# Jessica has asked for new dotplots of the 4C single cell data based on Sidu's
# interpretation of the ATAC data.
library(Seurat)
seu <- readRDS("~/Projects/sysAct_allFiles/DPayz_sysAct_seurat_20230327.rds")
levels(seu) <- c("Fibroblast I","Fibroblast II","Fibroblast III","Fibroblast IV","Fibroblast V","Fibroblast VI","Fibroblast VII","Fibroblast VIII","Epidermal I","Epidermal II","Epidermal III","Epidermal IV","Basal epidermis","Hematopoetic progenitor I","Hematopoetic progenitor II","Endothelial I","Endothelial II","RBC I","RBC II","Myeloid progenitor","Skeletal muscle" ,"T-cell")

genes <- as.data.frame(seu@assays$RNA@counts)
genes <- row.names(genes)

p53_set <- c(
  "TP53-AMEX60DD012903", # p53 (TP53)
  "WT1-AMEX60DDU001006343", # WT1
  "CLOCK-AMEX60DD045527", # CLOCK
  "CLOCK-AMEX60DD045528", # CLOCK
  "CLOCK-AMEX60DD045529", # CLOCK
  "PASD1-CLOCK-AMEX60DD037299", # CLOCK
  "LIM-1-LHX1-AMEX60DD054124", # LHX1
  "LHX3-AMEX60DD050656", # LHX3
  "POU2F3-AMEX60DD053808", # OCT11
  "POU5F1-AMEX60DD010762", # OCT11
  "POU2-POU5F1-AMEX60DD050901", # OCT11
  "POU2F1-AMEX60DD047046", # OCT11
  "POU2F2-AMEX60DD018055", # OCT11
  "POU3F1-AMEX60DD005829", # OCT6
  "SMAD2-AMEX60DD041720", # SMAD2
  "SMAD4-AMEX60DD015011", # SMAD4
  "SMAD4-AMEX60DD041763", # SMAD4
  "PAX5-AMEX60DD043936" #PAX5
  # OCT4 is the same loci as POU5F1-AMEX60DD010762 and POU2-POU5F1-AMEX60DD050901
)

DotPlot(seu,assay = "RNA",p53_set,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",p53_set,cols="RdBu") + theme(axis.text.x = element_text(angle = 90, hjust=1))

EN1_set <- c(
  "EN1-AMEX60DD056025", # En1
  "LHX9-AMEX60DD018728", # Lhx9
  "PITX1-AMEX60DD028957", # Pitx1
  "SMAD3-AMEX60DD003401", # Smad3
  "TWIST3-TWIST2-AMEX60DD029436", # Twist2
  "TWIST2-AMEX60DD055512" # Twist2
)
DotPlot(seu,assay = "RNA",EN1_set,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

AP1 <- c(
  # AP-1
  "C-FOS-FOS-AMEX60DD011151", #c-FOS/FOS
  "C-FOS-FOS-AMEX60DD011154", #c-FOS/FOS
  "FOSL1-AMEX60DD025233", # Fosl1
  "JUN-AMEX60DD019448", # Jun-AP1
  "JUN-AMEX60DD043789", # Jun-AP1
  "JUNB-AMEX60DD031925", # JunB
  "DLX2-AMEX60DD055671", # DLX2
  "DLX5-AMEX60DD022355", # DLX5
  "HOXA11-HOXA3-AMEX60DD021947", # Hoxa10
  "AB205-0109740-PAX6-AMEX60DD007993", # Pax6
  "PAX6-AMEX60DD004905", # Pax6
  "SIX2-AMEX60DD035940" # Six2
)

DotPlot(seu,assay = "RNA",sort(AP1, decreasing = T),cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

ATF <- c(
  "ATF1-AMEX60DD029673", # ATF-1 
  "ATF1-AMEX60DD035554", # ATF-1 
  "ATF2-AMEX60DD055620", # ATF-2
  "ATF3-AMEX60DDU001003554", # ATF-3
  "ATF4-AMEX60DD029230", # ATF-4
  "ATF5-AMEX60DD017394", # ATF-5
  "ATF6-AMEX60DD018411", # ATF-6
  "ATF6B-AMEX60DD010499", # ATF-6B
  "ATF7-AMEX60DD029531", # ATF-7
  "ATF7IP-AMEX60DD029096", # ATF7IP
  "ATF7IP2-AMEX60DD021479", # ATF7IP2
  "HEATR5B-AMEX60DD011156", # Possible BATF
  "HEATR5B-AMEX60DD011157", # Possible BATF
  "TTLL5-AMEX60DD011157", # Possible BATF
  "BATF2-AMEX60DD025272", # BATF-2
  "BATF3-AMEX60DD036067", # BATF-3
  "JDP2-AMEX60DD011155" # JDP2
)
DotPlot(seu,assay = "RNA",ATF,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

FOS <- c(
  # The AP-1 family 
  "C-FOS-FOS-AMEX60DD011151", #c-FOS
  "C-FOS-FOS-AMEX60DD011155", #c-FOS
  "FOSB-AMEX60DD024876", # FOSB
  "FOSB-AMEX60DD024877", # FOSB
  "FOSL1-AMEX60DD025233", # FRA-1
  "FOSL2-AMEX60DD036324" # FOSL2 FRA-2
)
DotPlot(seu,assay = "RNA",FOS,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

JUN <- c(
  # **JUN**
  "JUN-AMEX60DD019448", # Jun-AP1
  "JUN-AMEX60DD043789", # Jun-AP1
  "JUNB-AMEX60DD031925", # JunB
  "JUND-AMEX60DD014112" # JUND
)


DotPlot(seu,assay = "RNA",JUN,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

MAF <- c(
  # **MAF**
  "SIRT7-AMEX60DD031492", # MAFG? 
  "MAFF-MAFK-AMEX60DD029295",
  "MAFB-AMEX60DD009484",
  "MAFB-AMEX60DD028613"
)
DotPlot(seu,assay = "RNA",MAF,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
DLX <- c(
  # **Distal-Less Homeobox**
  "BOX5-DLX6-AMEX60DD009950", # DLX6
  "DLX6-AMEX60DD022357", # DLX6
  "DLX5-AMEX60DD022355", # DLX5
  "DLX3-AMEX60DD009952", # DLX3
  "DLX2-AMEX60DD055671", # DLX2
  "DLX1-AMEX60DD055670" # DLX1
)
DotPlot(seu,assay = "RNA",DLX,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

HOX <- c(  # **Homeobox**
  "CGPICR-010019-SHOX2-AMEX60DD001828",
  "SHOX-AMEX60DD047789",
  "PHOX2B-AMEX60DD045633",
  "HOXD13-AMEX60DD055612",
  "HOXD11-AMEX60DD055611",
  "HOXD10-AMEX60DD055610",
  "HOXD4-AMEX60DD055607",
  "HOXD3-AMEX60DD055607",
  "AB205-0035640-HOXD1-AMEX60DD055606",
  "HOXC13-AMEX60DD029508",
  "HOXC11-AMEX60DD029495",
  "D130011F21RIK-HOXC8-AMEX60DD029496",
  "HOXC5-AMEX60DD029497",
  "HOXC4-AMEX60DD029499",
  "HOXC3-HOXD3-AMEX60DD029486",
  "HOXB13-AMEX60DD010173",
  "HOXB7-AMEX60DD010179",
  "HOXB3-AMEX60DD010179",
  "HOXB2-AMEX60DD010181",
  "AB205-0166220-HOXB1-AMEX60DD010182",
  "HOXA13-AMEX60DD021953",
  "HOXA11-HOXA3-AMEX60DD021947",
  "HOXA2-AMEX60DD021945"
)
DotPlot(seu,assay = "RNA",HOX,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

PAX <- c(
  # **PAX**
  "PAX8-AMEX60DD001281",
  "PAX6-AMEX60DD004905",
  "PAX9-AMEX60DD012252",
  "PAXIP1-AMEX60DD021590",
  "PAXIP1-AMEX60DD021592",
  "PAX1-AMEX60DD035469",
  "PAX5-AMEX60DD043936",
  "DV515-00017576-PAXBP1-AMEX60DD047311",
  "PAX2-AMEX60DD051628",
  "PAXX-AMEX60DD050906",
  "PAX7-AMEX60DD052322",
  "AB205-0109740-PAX6-AMEX60DD007993"
)
DotPlot(seu,assay = "RNA",PAX,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
SIX <- c(
  #**SIX**
  "SIX5-AMEX60DD024752",
  "SIX4-AMEX60DD011855",
  "SIX2-AMEX60DD035940",
  "SIX1-AMEX60DD024753",
  "SIX1-AMEX60DD011856"
)
DotPlot(seu,assay = "RNA",SIX,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
MAF <- c(
  # **MAF**
  "SIRT7-AMEX60DD031492", # MAFG? 
  "MAFF-MAFK-AMEX60DD029295",
  "MAFB-AMEX60DD009484",
  "MAFB-AMEX60DD028613",
  "MAF1-AMEX60DD040812"
)
DotPlot(seu,assay = "RNA",MAF,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

nMYC <- c(
  "NDRG1-AMEX60DD040549", # n-Myc
  "NDRG1-AMEX60DD044627", # n-Myc
  "ZNF683-PRDM1-AMEX60DD005654", # Prdm1
  "PRDM1-AMEX60DD034088", # Prdm1
  "POU2F2-AMEX60DD018055", # Oct2 POU2F2
  "ATF3-AMEX60DDU001003554" # ATF-3
) 
DotPlot(seu,assay = "RNA",sort(nMYC, decreasing = T),cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

NDRG <- c(
# NDRG family
"NDRG4-AMEX60DD015507", # n-Myc
"NDRG3-AMEX60DD028553", # n-Myc
"NDRG2-AMEX60DD008501", # n-Myc
"NDRG1-AMEX60DD040549", # n-Myc
"NDRG1-AMEX60DD044627" # n-Myc
)
DotPlot(seu,assay = "RNA",NDRG,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

PRDM <- c(
# PRDM family
"PRDM16-AMEX60DD051365",
"PRDM15-AMEX60DD047134",
"PRDM12-AMEX60DD050368",
"PRDM11-AMEX60DD005203",
"EOD39-8286-PRDM11-AMEX60DD040334",
"PRDM10-AMEX60DD053575",
"ZNF683-PRDM1-AMEX60DD005654",
"PRDM1-AMEX60DD034088"
)
DotPlot(seu,assay = "RNA",PRDM,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

HOXnSIX <- c(
  "HOXD13-AMEX60DD055612",
  # Hoxd12 Not present
  "HOXB13-AMEX60DD010173",
  "HOXA13-AMEX60DD021953",
  "HOXA11-HOXA3-AMEX60DD021947",
  "MEIS1-AMEX60DD035598", # Meis1
  "SIX1-AMEX60DD011856", # Six1
  "SIX1-AMEX60DD024753" # Six1
)
DotPlot(seu,assay = "RNA",sort(HOXnSIX, decreasing = T),cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

LEPTIN <- c(
  "LEP-AMEX60DD007977", # Leptin
  "JUN-AMEX60DD019448", # c-Jun
  "JUN-AMEX60DD043789", # c-Jun
  "C-MYC-MYC-AMEX60DD040490" #c-Myc
)
DotPlot(seu,assay = "RNA",sort(LEPTIN, decreasing = T),cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

patterns <- c("AMEX60DD000569",
             "AMEX60DD000570",
             "AMEX60DD000963",
             "AMEX60DD001255",
             "AMEX60DD001368",
             "AMEX60DD001369",
             "AMEX60DD001640",
             "AMEX60DD002028",
             "AMEX60DD002376",
             "AMEX60DD002437",
             "AMEX60DD002521",
             "AMEX60DD002668",
             "AMEX60DD002708",
             "AMEX60DD002709",
             "AMEX60DD002999",
             "AMEX60DD003650",
             "AMEX60DD003651",
             "AMEX60DD004113",
             "AMEX60DD004724",
             "AMEX60DD004725",
             "AMEX60DD004753",
             "AMEX60DD004771",
             "AMEX60DD004772",
             "AMEX60DD004832",
             "AMEX60DD004833",
             "AMEX60DD005012",
             "AMEX60DD005013",
             "AMEX60DD005142",
             "AMEX60DD005314",
             "AMEX60DD005411",
             "AMEX60DD005422",
             "AMEX60DD005525",
             "AMEX60DD005547",
             "AMEX60DD005550",
             "AMEX60DD005551",
             "AMEX60DD005698",
             "AMEX60DD005819",
             "AMEX60DD005942",
             "AMEX60DD006052",
             "AMEX60DD006178",
             "AMEX60DD006179",
             "AMEX60DD006290",
             "AMEX60DD006292",
             "AMEX60DD006373",
             "AMEX60DD006480",
             "AMEX60DD006538",
             "AMEX60DD006750",
             "AMEX60DD006751",
             "AMEX60DD006761",
             "AMEX60DD006799",
             "AMEX60DD006800",
             "AMEX60DD006806",
             "AMEX60DD006807",
             "AMEX60DD006808",
             "AMEX60DD006931",
             "AMEX60DD006976",
             "AMEX60DD006977",
             "AMEX60DD006998",
             "AMEX60DD006999",
             "AMEX60DD007000",
             "AMEX60DD007039",
             "AMEX60DD007040",
             "AMEX60DD007595",
             "AMEX60DD007596",
             "AMEX60DD007740",
             "AMEX60DD008148",
             "AMEX60DD008764",
             "AMEX60DD008788",
             "AMEX60DD008850",
             "AMEX60DD009090",
             "AMEX60DD009153",
             "AMEX60DD009155",
             "AMEX60DD009176",
             "AMEX60DD009181",
             "AMEX60DD009445",
             "AMEX60DD009484",
             "AMEX60DD009721",
             "AMEX60DD009723",
             "AMEX60DD010057",
             "AMEX60DD010206",
             "AMEX60DD010209",
             "AMEX60DD010287",
             "AMEX60DD010764",
             "AMEX60DD010836",
             "AMEX60DD010837",
             "AMEX60DD010971",
             "AMEX60DD011028",
             "AMEX60DD011029",
             "AMEX60DD011139",
             "AMEX60DD011302",
             "AMEX60DD011326",
             "AMEX60DD011327",
             "AMEX60DD011696",
             "AMEX60DD011699",
             "AMEX60DD011915",
             "AMEX60DD012004",
             "AMEX60DD012005",
             "AMEX60DD012170",
             "AMEX60DD012175",
             "AMEX60DD012219",
             "AMEX60DD012229",
             "AMEX60DD012261",
             "AMEX60DD012276",
             "AMEX60DD012772",
             "AMEX60DD012880",
             "AMEX60DD013102",
             "AMEX60DD013328",
             "AMEX60DD013331",
             "AMEX60DD013332",
             "AMEX60DD013929",
             "AMEX60DD014112",
             "AMEX60DD014386",
             "AMEX60DD014528",
             "AMEX60DD014529",
             "AMEX60DD014530",
             "AMEX60DD014531",
             "AMEX60DD014532",
             "AMEX60DD014533",
             "AMEX60DD014944",
             "AMEX60DD015023",
             "AMEX60DD015258",
             "AMEX60DD015450",
             "AMEX60DD015649",
             "AMEX60DD015650",
             "AMEX60DD015651",
             "AMEX60DD015664",
             "AMEX60DD015765",
             "AMEX60DD015832",
             "AMEX60DD015833",
             "AMEX60DD015835",
             "AMEX60DD015946",
             "AMEX60DD016289",
             "AMEX60DD017050",
             "AMEX60DD017233",
             "AMEX60DD017266",
             "AMEX60DD017302",
             "AMEX60DD017592",
             "AMEX60DD017651",
             "AMEX60DD017652",
             "AMEX60DD017653",
             "AMEX60DD017654",
             "AMEX60DD017655",
             "AMEX60DD017656",
             "AMEX60DD017803",
             "AMEX60DD017820",
             "AMEX60DD018277",
             "AMEX60DD018549",
             "AMEX60DD018550",
             "AMEX60DD019077",
             "AMEX60DD019394",
             "AMEX60DD019448",
             "AMEX60DD019552",
             "AMEX60DD019553",
             "AMEX60DD019621",
             "AMEX60DD020289",
             "AMEX60DD020295",
             "AMEX60DD020296",
             "AMEX60DD020747",
             "AMEX60DD020875",
             "AMEX60DD020982",
             "AMEX60DD020987",
             "AMEX60DD020988",
             "AMEX60DD021233",
             "AMEX60DD021382",
             "AMEX60DD021993",
             "AMEX60DD022215",
             "AMEX60DD022766",
             "AMEX60DD022767",
             "AMEX60DD022768",
             "AMEX60DD022769",
             "AMEX60DD022770",
             "AMEX60DD023136",
             "AMEX60DD023137",
             "AMEX60DD023229",
             "AMEX60DD023508",
             "AMEX60DD023837",
             "AMEX60DD023909",
             "AMEX60DD023939",
             "AMEX60DD023999",
             "AMEX60DD024365",
             "AMEX60DD024366",
             "AMEX60DD024704",
             "AMEX60DD024803",
             "AMEX60DD024834",
             "AMEX60DD025124",
             "AMEX60DD025372",
             "AMEX60DD025398",
             "AMEX60DD025428",
             "AMEX60DD025660",
             "AMEX60DD025661",
             "AMEX60DD025953",
             "AMEX60DD025954",
             "AMEX60DD025955",
             "AMEX60DD026039",
             "AMEX60DD026223",
             "AMEX60DD026763",
             "AMEX60DD026850",
             "AMEX60DD026909",
             "AMEX60DD026958",
             "AMEX60DD026959",
             "AMEX60DD026960",
             "AMEX60DD026961",
             "AMEX60DD027046",
             "AMEX60DD027047",
             "AMEX60DD027048",
             "AMEX60DD027152",
             "AMEX60DD027155",
             "AMEX60DD027344",
             "AMEX60DD027504",
             "AMEX60DD027592",
             "AMEX60DD027593",
             "AMEX60DD027594",
             "AMEX60DD027595",
             "AMEX60DD027906",
             "AMEX60DD027908",
             "AMEX60DD027909",
             "AMEX60DD028015",
             "AMEX60DD028016",
             "AMEX60DD028017",
             "AMEX60DD028018",
             "AMEX60DD028458",
             "AMEX60DD028613",
             "AMEX60DD028770",
             "AMEX60DD028785",
             "AMEX60DD028786",
             "AMEX60DD028921",
             "AMEX60DD029051",
             "AMEX60DD029172",
             "AMEX60DD029295",
             "AMEX60DD029302",
             "AMEX60DD029478",
             "AMEX60DD029551",
             "AMEX60DD029552",
             "AMEX60DD029553",
             "AMEX60DD029558",
             "AMEX60DD029622",
             "AMEX60DD029700",
             "AMEX60DD029784",
             "AMEX60DD029846",
             "AMEX60DD029859",
             "AMEX60DD029860",
             "AMEX60DD030034",
             "AMEX60DD030037",
             "AMEX60DD030070",
             "AMEX60DD030440",
             "AMEX60DD030500",
             "AMEX60DD030501",
             "AMEX60DD030603",
             "AMEX60DD030604",
             "AMEX60DD030605",
             "AMEX60DD030832",
             "AMEX60DD031066",
             "AMEX60DD031117",
             "AMEX60DD031118",
             "AMEX60DD031119",
             "AMEX60DD031120",
             "AMEX60DD031121",
             "AMEX60DD031160",
             "AMEX60DD031236",
             "AMEX60DD031299",
             "AMEX60DD031492",
             "AMEX60DD031538",
             "AMEX60DD031634",
             "AMEX60DD031642",
             "AMEX60DD031647",
             "AMEX60DD031821",
             "AMEX60DD031925",
             "AMEX60DD031989",
             "AMEX60DD032013",
             "AMEX60DD032172",
             "AMEX60DD032637",
             "AMEX60DD032807",
             "AMEX60DD032826",
             "AMEX60DD032827",
             "AMEX60DD032871",
             "AMEX60DD032918",
             "AMEX60DD032920",
             "AMEX60DD032921",
             "AMEX60DD033220",
             "AMEX60DD033221",
             "AMEX60DD033222",
             "AMEX60DD033320",
             "AMEX60DD033321",
             "AMEX60DD033458",
             "AMEX60DD033459",
             "AMEX60DD033523",
             "AMEX60DD033563",
             "AMEX60DD033941",
             "AMEX60DD033944",
             "AMEX60DD034172",
             "AMEX60DD034349",
             "AMEX60DD034390",
             "AMEX60DD034391",
             "AMEX60DD034504",
             "AMEX60DD034577",
             "AMEX60DD034578",
             "AMEX60DD034656",
             "AMEX60DD034772",
             "AMEX60DD034773",
             "AMEX60DD034864",
             "AMEX60DD035450",
             "AMEX60DD035760",
             "AMEX60DD035892",
             "AMEX60DD035893",
             "AMEX60DD035896",
             "AMEX60DD035897",
             "AMEX60DD035936",
             "AMEX60DD036019",
             "AMEX60DD036329",
             "AMEX60DD036784",
             "AMEX60DD036816",
             "AMEX60DD036817",
             "AMEX60DD036841",
             "AMEX60DD037080",
             "AMEX60DD037081",
             "AMEX60DD037082",
             "AMEX60DD037083",
             "AMEX60DD037084",
             "AMEX60DD037085",
             "AMEX60DD037199",
             "AMEX60DD037200",
             "AMEX60DD037201",
             "AMEX60DD037226",
             "AMEX60DD037227",
             "AMEX60DD037395",
             "AMEX60DD037450",
             "AMEX60DD037561",
             "AMEX60DD037749",
             "AMEX60DD037864",
             "AMEX60DD038107",
             "AMEX60DD038319",
             "AMEX60DD038379",
             "AMEX60DD038381",
             "AMEX60DD038386",
             "AMEX60DD038468",
             "AMEX60DD038477",
             "AMEX60DD038478",
             "AMEX60DD038738",
             "AMEX60DD038834",
             "AMEX60DD038928",
             "AMEX60DD039000",
             "AMEX60DD039073",
             "AMEX60DD039435",
             "AMEX60DD039494",
             "AMEX60DD039561",
             "AMEX60DD039712",
             "AMEX60DD039713",
             "AMEX60DD039720",
             "AMEX60DD039791",
             "AMEX60DD039921",
             "AMEX60DD039922",
             "AMEX60DD040024",
             "AMEX60DD040395",
             "AMEX60DD040396",
             "AMEX60DD040397",
             "AMEX60DD040678",
             "AMEX60DD040733",
             "AMEX60DD041668",
             "AMEX60DD041795",
             "AMEX60DD041796",
             "AMEX60DD041852",
             "AMEX60DD041908",
             "AMEX60DD042370",
             "AMEX60DD042419",
             "AMEX60DD042420",
             "AMEX60DD042452",
             "AMEX60DD042453",
             "AMEX60DD042466",
             "AMEX60DD042467",
             "AMEX60DD043101",
             "AMEX60DD043102",
             "AMEX60DD043103",
             "AMEX60DD043128",
             "AMEX60DD043789",
             "AMEX60DD044800",
             "AMEX60DD044801",
             "AMEX60DD044802",
             "AMEX60DD045262",
             "AMEX60DD045265",
             "AMEX60DD045294",
             "AMEX60DD045929",
             "AMEX60DD046035",
             "AMEX60DD046138",
             "AMEX60DD046404",
             "AMEX60DD046676",
             "AMEX60DD046871",
             "AMEX60DD046872",
             "AMEX60DD046924",
             "AMEX60DD046925",
             "AMEX60DD047184",
             "AMEX60DD047360",
             "AMEX60DD047530",
             "AMEX60DD047838",
             "AMEX60DD048270",
             "AMEX60DD048281",
             "AMEX60DD048355",
             "AMEX60DD048489",
             "AMEX60DD048490",
             "AMEX60DD048491",
             "AMEX60DD048493",
             "AMEX60DD048494",
             "AMEX60DD048504",
             "AMEX60DD048613",
             "AMEX60DD048974",
             "AMEX60DD049073",
             "AMEX60DD049075",
             "AMEX60DD049077",
             "AMEX60DD049150",
             "AMEX60DD049295",
             "AMEX60DD049895",
             "AMEX60DD050428",
             "AMEX60DD050429",
             "AMEX60DD050434",
             "AMEX60DD050724",
             "AMEX60DD051347",
             "AMEX60DD051348",
             "AMEX60DD051351",
             "AMEX60DD051448",
             "AMEX60DD051458",
             "AMEX60DD051577",
             "AMEX60DD051579",
             "AMEX60DD051582",
             "AMEX60DD052489",
             "AMEX60DD052490",
             "AMEX60DD052551",
             "AMEX60DD052712",
             "AMEX60DD052723",
             "AMEX60DD052740",
             "AMEX60DD052855",
             "AMEX60DD052913",
             "AMEX60DD052914",
             "AMEX60DD053002",
             "AMEX60DD053004",
             "AMEX60DD053005",
             "AMEX60DD053152",
             "AMEX60DD053743",
             "AMEX60DD053863",
             "AMEX60DD054110",
             "AMEX60DD054236",
             "AMEX60DD054285",
             "AMEX60DD054860",
             "AMEX60DD054873",
             "AMEX60DD054874",
             "AMEX60DD055651",
             "AMEX60DD055695",
             "AMEX60DD055696",
             "AMEX60DD055701",
             "AMEX60DD055801",
             "AMEX60DD056120",
             "AMEX60DDU001000911",
             "AMEX60DDU001001824",
             "AMEX60DDU001001832",
             "AMEX60DDU001001834",
             "AMEX60DDU001001835",
             "AMEX60DDU001001836",
             "AMEX60DDU001001839",
             "AMEX60DDU001001841",
             "AMEX60DDU001001842",
             "AMEX60DDU001001846",
             "AMEX60DDU001001857",
             "AMEX60DDU001001858",
             "AMEX60DDU001001859",
             "AMEX60DDU001002453",
             "AMEX60DDU001002454",
             "AMEX60DDU00100576",
             "AMEX60DDU001005765",
             "AMEX60DDU001005766",
             "AMEX60DDU001005767",
             "AMEX60DDU001005768",
             "AMEX60DDU001005769",
             "AMEX60DDU001006136",
             "AMEX60DDU001006586",
             "AMEX60DDU001008595",
             "AMEX60DDU001008614",
             "AMEX60DDU001009496",
             "AMEX60DDU001009497",
             "AMEX60DDU001010229",
             "AMEX60DDU001010327",
             "AMEX60DDU001010331",
             "AMEX60DDU001010332",
             "AMEX60DDU001010333",
             "AMEX60DDU001010334",
             "AMEX60DDU001010828",
             "AMEX60DDU001012929",
             "AMEX60DDU001012930",
             "AMEX60DDU001014364",
             "AMEX60DDU00101810",
             "AMEX60DDU001018100",
             "AMEX60DDU001018101",
             "AMEX60DDU001018102",
             "AMEX60DDU001018103",
             "AMEX60DDU001018104",
             "AMEX60DDU001018105",
             "AMEX60DDU001023251",
             "AMEX60DDU001025731",
             "AMEX60DDU001030413",
             "AMEX60DDU001035578",
             "AMEX60DDU00103579",
             "AMEX60DDU001035791",
             "AMEX60DDU001035792",
             "AMEX60DDU001035793",
             "AMEX60DDU001035794",
             "AMEX60DDU001035796",
             "AMEX60DDU001036217",
             "AMEX60DDU001037208",
             "AMEX60DDU001037209",
             "AMEX60DDU001037210",
             "AMEX60DDU001037211",
             "AMEX60DDU001038740")
# foo <- foo %>% distinct()
txf_rows <- list()
for (pattern in patterns) {
filtered_rows <- genes[grepl(pattern, genes$genes), ]  # Using grepl for example
txf_rows[[length(txf_rows) + 1]] <- filtered_rows
}
txf_genes <- matrix(unlist(txf_rows), nrow = length(txf_rows), byrow = TRUE)
txf_genes <- sort(txf_genes, decreasing = T)
txf_genes <- as.data.frame(txf_genes)
txf_genes <- txf_genes %>% distinct()
foo <- as.list(txf_genes)
DotPlot(seu,assay = "RNA",sort(txf_genes$txf_genes[451:467], decreasing = T),cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[31:60],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[61:90],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[91:120],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[121:150],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[151:180],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[181:210],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[211:240],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[241:270],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[271:300],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[301:330],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[331:360],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[361:390],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[391:420],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[421:450],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# DotPlot(seu,assay = "RNA",txf_genes[451:467],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()
# # Pattern to find (example)
pattern <- "-"

# Invert the selection to get rows without "-"

known <- txf_genes[ grepl(pattern, txf_genes$txf_genes), ]
known <- sort(known, decreasing = T)

unknown <- txf_genes[ !grepl(pattern, txf_genes$txf_genes), ]
unknown <- sort(unknown, decreasing = T)

DotPlot(seu,assay = "RNA",known[1:100],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()

DotPlot(seu,assay = "RNA",unknown[31:60],cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()




## Create PCA plot?

# rm(list = ls())

library("ChIPseeker")
# library("clusterProfiler")
library(GenomicFeatures)
library(GenomicRanges)
library(stats)
# library(gridExtra)
# library(ATACseqQC)
library(dplyr)
library(DESeq2)
# library(BSgenome.Amexicanum.NCBI.ambMex60DD)
library(ggplot2)
library(ggrepel)
library(tidyverse)

setwd("~/Projects/sysAct_allFiles/ATAC_Mar5_2024")
txdb <-
  makeTxDbFromGFF(file = "/Users/sjblair/Projects/sysAct_allFiles/cluster_sysAct_files/4NSAC_UPDATE20231122/amexG.chunked.annot.sorted.gtf", format = "gtf")

# counts <- "/Users/sjblair/Projects/sysAct_allFiles/ATAC_Mar5_2024/ATAC_COUNTS_Mar13_2024"


counts <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/bl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","bl1"))
counts[1:5] <- NULL
bl4 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/bl4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","bl4"))

cl4 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/cl4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","cl4"))
cl5 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/cl5.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","cl5"))

in3 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/in3.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","in3"))
in4 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/in4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","in4"))
in5 <- read.table("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/in5.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","in5"))

counts$bl4 <- bl4$bl4
counts$cl4 <- cl4$cl4
counts$cl5 <- cl5$cl5
counts$in3 <- in3$in3
counts$in4 <- in4$in4
counts$in5 <- in5$in5
rm(bl2,bl3,bl4,in1,in2,in3,in4,in5,cl3,cl4,cl5)


pca <- prcomp(t(counts), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]

# Create the plot
ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(counts)), size = 5) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, scaled and centered") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label =colnames(counts)), size = 5,nudge_x = 5,nudge_y = 1) +
  theme_classic()
cbinder <- read.table("/Users/sjblair/Projects/sysAct_allFiles/ATAC_Mar5_2024/ATAC_COUNTS_Mar15_2024/bl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","bl1"))
cbinder[4:6] <- NULL
cbinder$peak <- rownames(cbinder)
cbinder$grange <- paste0(cbinder$Chr,":",cbinder$Start,"-",cbinder$End)

anno <- read.table("/Users/sjblair/Projects/sysAct_allFiles/ATAC_Mar5_2024/ATAC_COUNTS_Mar15_2024/bl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl1"))
anno[4:6] <- NULL
anno$peak <- rownames(anno)
anno <- makeGRangesFromDataFrame(anno)

anno <- annotatePeak(anno,tssRegion=c(-3000, 3000),TxDb=txdb) 
anno <- as.data.frame(anno)
anno$grange <- paste0(anno$seqnames,":",anno$start,"-",anno$end)

bindanno <- left_join(anno, cbinder)
rownames <- bindanno
rownames[1:12] <- NULL
rownames[2:6] <- NULL

rownames$transcriptId <- gsub(" \\[","..",rownames$transcriptId)
rownames$transcriptId <- gsub("\\]","..",rownames$transcriptId)
rownames$transcriptId <- gsub("..nr..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("..hs..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("..hs..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("\\|","^",rownames$transcriptId)
rownames$transcriptId <- gsub("\\^\\^","^",rownames$transcriptId)
rownames$rownames <- paste0(rownames$transcriptId,"^",rownames$peak)

rownames$transcriptId <- NULL
counts$peak <- rownames(counts)

cts <- left_join(rownames, counts)
cts$peak <- NULL
rownames(cts) <- cts$rownames
cts$rownames <- NULL

in_cl.cts <- cts[4:11]


# compare only the conditions from before
foo <- cts
foo$gene <- row.names(foo)
foo[1:4] <- NULL
foo[5:8] <- NULL

coldata <- colnames(foo)
condition <- c("Activated","Activated","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
foo.dds <- DESeqDataSetFromMatrix(foo, colData = coldata,
                                  design = ~ condition)
foo.dds <- DESeq(foo.dds)
foo.res <- results(foo.dds, contrast=c("condition","Activated","Homeostatic"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

foo.res <- foo.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

foo.pca <- prcomp(t(foo), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- foo.pca$x[, 1]
pc2 <- foo.pca$x[, 2]

# Create the plot
in12_cl45_pcaPlot <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(foo)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, blastema omit, scaled and centered") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(foo)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()


# Brian said on Friday to compare cl3,4,5 vs intact all and then 
# Then for comparison run the same with dropping in4 to see how much that changes results
# cl3,4,5 vs intact all has 0 results, I showed above that in1,2 vs cl4,5 has 50.  
# Let's see whaet removing in4 from the matrix.
# for comparision, in_cl results for FOXO1 are as follows.
in_cl.res[39126, ]
# FOXO1^AMEX60DD301049150.1^Peak_39151     245.          0.198 0.133  1.49  0.137 0.839
in1235_cl345.cts <- in_cl.cts
in1235_cl345.cts$in4 <- NULL

coldata <- colnames(in1235_cl345.cts)
condition <- c("Activated","Activated","Activated","Homeostatic","Homeostatic","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
in1235_cl345.dds <- DESeqDataSetFromMatrix(in1235_cl345.cts, colData = coldata,
                                           design = ~ condition)
in1235_cl345.dds <- DESeq(in1235_cl345.dds)
in1235_cl345.res <- results(in1235_cl345.dds, contrast=c("condition","Activated","Homeostatic"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

in1235_cl345.res <- in1235_cl345.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

in1235_cl345.pca <- prcomp(t(in1235_cl345.cts), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- in1235_cl345.pca$x[, 1]
pc2 <- in1235_cl345.pca$x[, 2]

# Create the plot
in1235_cl345_pcaPlot <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(in1235_cl345.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, in1,2,3,5 vs cl 3,4,5, scaled and centered") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(in1235_cl345.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()

## plot the TSS trend
# ggplot(anno[14], aes(factor(anno[14],max(anno$distanceToTSS) / min(anno$distanceToTSS), group = 1)) +
#   geom_point() +
#   geom_line() +
#   xlab("Year published") +
#   ylab("Proportion of OA full-texts in Europe PMC")

save.image(file = "ATAC_Mar182024.rdata")

in1235_cl345.loading <- in1235_cl345.pca$rotation
in1235_cl345.loading[39126, ]
in1235_cl345.cts[39126, ]

in1235_cl345.loading[head(order(in1235_cl345.loading$PC2, decreasing = T),n=20),]
in1235_cl345.loading[head(order(in1235_cl345.loading$PC2, decreasing = F),n=20),]

in1235_cl345.loading[head(order(in1235_cl345.loading$PC2, decreasing = T)),]
in1235_cl345.loading[head(order(in1235_cl345.loading$PC2, decreasing = F)),]


# next I will create a subset of the counts that are only at TSS distance of 0.

in1235_cl345.TSS0 <- bindanno %>%
  data.frame()  %>% 
  filter(distanceToTSS == 0)

in1235_cl345.TSS0$transcriptId <- gsub(" \\[","..",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("\\]","..",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("..nr..","^",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("..hs..","^",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("..hs..","^",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("\\|","^",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$transcriptId <- gsub("\\^\\^","^",in1235_cl345.TSS0$transcriptId)
in1235_cl345.TSS0$rownames <- paste0(in1235_cl345.TSS0$transcriptId,"^",in1235_cl345.TSS0$peak)

in1235_cl345.TSS0[1:18] <- NULL


in1235_cl345.TSS0 <- left_join(in1235_cl345.TSS0, counts)
rownames(in1235_cl345.TSS0) <- in1235_cl345.TSS0$rownames
in1235_cl345.TSS0[1:2] <- NULL
in1235_cl345.TSS0[1:3] <- NULL
in1235_cl345.TSS0$in4 <- NULL

save.image(file = "ATAC_Mar182024.rdata")
in1235_cl345.TSS0.cts <- in1235_cl345.TSS0
coldata <- colnames(in1235_cl345.TSS0.cts)
condition <- c("Activated","Activated","Activated","Homeostatic","Homeostatic","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
in1235_cl345.TSS0.dds <- DESeqDataSetFromMatrix(in1235_cl345.TSS0.cts, colData = coldata,
                                                design = ~ condition)
in1235_cl345.TSS0.dds <- DESeq(in1235_cl345.TSS0.dds)
in1235_cl345.TSS0.res <- results(in1235_cl345.TSS0.dds, contrast=c("condition","Activated","Homeostatic"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

in1235_cl345.TSS0.res <- in1235_cl345.TSS0.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene")
filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

in1235_cl345.TSS0.pca <- prcomp(t(in1235_cl345.TSS0.cts), center = T, scale. = T)

# Select the first two principal components
pc1 <- in1235_cl345.TSS0.pca$x[, 1]
pc2 <- in1235_cl345.TSS0.pca$x[, 2]

# Create the plot
in1235_cl345.TSS0_pcaPlot <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(in1235_cl345.TSS0.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, in1,2,3,5 vs cl 3,5, TSS Adjacent") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(in1235_cl345.TSS0.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()


# let's see what happens when we remove intact 5
in123_cl345.TSS0.cts <- in1235_cl345.TSS0.cts
in123_cl345.TSS0.cts$in5 <- NULL
coldata <- colnames(in123_cl345.TSS0.cts)
condition <- c("Activated","Activated","Activated","Homeostatic","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
in123_cl345.TSS0.dds <- DESeqDataSetFromMatrix(in123_cl345.TSS0.cts, colData = coldata,
                                               design = ~ condition)
in123_cl345.TSS0.dds <- DESeq(in123_cl345.TSS0.dds)
in123_cl345.TSS0.res <- results(in123_cl345.TSS0.dds, contrast=c("condition","Activated","Homeostatic"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

in123_cl345.TSS0.res <- in123_cl345.TSS0.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

in123_cl345.TSS0.pca <- prcomp(t(in123_cl345.TSS0.cts), center = T, scale. = T)

# Select the first two principal components
pc1 <- in123_cl345.TSS0.pca$x[, 1]
pc2 <- in123_cl345.TSS0.pca$x[, 2]

# Create the plot
in123_cl345.TSS0_pcaPlot <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(in123_cl345.TSS0.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, in1,2,3,5 vs cl 3,5, TSS Adjacent") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(in123_cl345.TSS0.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()

# cl vs bl
bl123_cl345.cts <- cts
bl123_cl345.cts[7:11] <- NULL
coldata <- colnames(bl123_cl345.cts)
condition <- c("Blastema","Blastema","Blastema","Activated","Activated","Activated")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
bl123_cl345.dds <- DESeqDataSetFromMatrix(bl123_cl345.cts, colData = coldata,
                                          design = ~ condition)
bl123_cl345.dds <- DESeq(bl123_cl345.dds)
bl123_cl345.res <- results(bl123_cl345.dds, contrast=c("condition","Activated","Blastema"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

bl123_cl345.res <- bl123_cl345.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

bl123_cl345.pca <- prcomp(t(bl123_cl345.cts), center = T, scale. = T)

# Select the first two principal components
pc1 <- bl123_cl345.pca$x[, 1]
pc2 <- bl123_cl345.pca$x[, 2]

# Create the plot
bl123_cl345.pca <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(bl123_cl345.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, bl123 vs cl345") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(bl123_cl345.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()

write.table(bl123_cl345.res, file = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/PCA_bl123_cl345.res.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

# blastema vs intact
bl123_in1235.cts <- cts
bl123_in1235.cts[4:6] <- NULL
bl123_in1235.cts[7] <- NULL

coldata <- colnames(bl123_in1235.cts)
condition <- c("Blastema","Blastema","Blastema","Homeostatic","Homeostatic","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
bl123_in1235.dds <- DESeqDataSetFromMatrix(bl123_in1235.cts, colData = coldata,
                                           design = ~ condition)
bl123_in1235.dds <- DESeq(bl123_in1235.dds)
bl123_in1235.res <- results(bl123_in1235.dds, contrast=c("condition","Homeostatic","Blastema"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

bl123_in1235.res <- bl123_in1235.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
# that is interesting, 50 results between these conditions used from last year.  

bl123_in1235.pca <- prcomp(t(bl123_in1235.cts), center = T, scale. = T)

# Select the first two principal components
pc1 <- bl123_in1235.pca$x[, 1]
pc2 <- bl123_in1235.pca$x[, 2]

# Create the plot
bl123_in1235.pca <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(bl123_in1235.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC counts, bl123 vs in1235") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(bl123_in1235.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()

write.table(bl123_in1235.res, file = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/PCA_bl123_in1235.res.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

# volcano plot
library(EnhancedVolcano)


EnhancedVolcano(bl123_cl345.res,
                lab = gsub("\\^AMEX.*","",bl123_cl345.res$gene),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('VCAM1','KCTD12','ADAM12',
                #              'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                title = 'ATAC results, Activated vs Blastema',
                subtitle = "",
                labSize = 2,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = T,
                colAlpha = 4/5,
                legendPosition = 'none',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = T,
                widthConnectors = 0.5,
                colConnectors = 'black',
                max.overlaps = 20)

EnhancedVolcano(bl123_in1235.res,
                lab = gsub("\\^AMEX.*","",bl123_in1235.res$gene),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('VCAM1','KCTD12','ADAM12',
                #              'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                title = 'ATAC results, Homeostatic vs Blastema',
                subtitle = "",
                labSize = 2,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = T,
                colAlpha = 4/5,
                legendPosition = 'none',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = T,
                widthConnectors = 0.5,
                colConnectors = 'black',
                max.overlaps = 20)

# Lets do a pathway on these results
bl123_cl345.pathway <- bl123_cl345.res[c("gene","log2FoldChange","padj")]
bl123_cl345.pathway$gene <- gsub("\\^.*","",bl123_cl345.pathway$gene)
library(pathfindR)

bl123_cl345_pathway <- run_pathfindR(
  input = bl123_cl345.pathway,
  convert2alias = FALSE,
  gene_sets = "KEGG",
  # custom_genes = kegg_genes,
  # custom_descriptions = kegg_descriptions,
  iterations = 25,
  n_processes = 8,
  output_dir = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_cl_bl"
) # Found 580/1647 genes, found 107 KEGG pathways
enrichment_chart(result_df = bl123_cl345_pathway, top_terms = 107)

write.table(bl123_cl345_pathway,"/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_cl_bl_Enrichment_KEGG.tsv",quote = F, sep = '\t',row.names = F, col.names = T)

bl123_cl345_pathway_GO <- run_pathfindR(
  input = bl123_cl345.pathway,
  convert2alias = FALSE,
  gene_sets = "GO-All",
  # custom_genes = kegg_genes,
  # custom_descriptions = kegg_descriptions,
  iterations = 25,
  n_processes = 8,
  output_dir = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_cl_bl_GO"
) # Found 580/1647 genes, found 123 GO-All terms
write.table(bl123_cl345_pathway,"/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_cl_bl_Enrichment_GO.tsv",quote = F, sep = '\t',row.names = F, col.names = T)

enrichment_chart(result_df = bl123_cl345_pathway_GO, top_terms = 123)

# now to do the same with intact against blastema
# Lets do a pathway on these results
bl123_in1235.pathway <- bl123_in1235.res[c("gene","log2FoldChange","padj")]
bl123_in1235.pathway$gene <- gsub("\\^.*","",bl123_in1235.pathway$gene)

bl123_in1235_pathway_KEGG <- run_pathfindR(
  input = bl123_in1235.pathway,
  convert2alias = FALSE,
  gene_sets = "KEGG",
  # custom_genes = kegg_genes,
  # custom_descriptions = kegg_descriptions,
  iterations = 25,
  n_processes = 8,
  output_dir = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_in_bl_KEGG"
) # Found 372/1036 genes, found 78 KEGG pathways

write.table(bl123_in1235_pathway_KEGG,"/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_in_bl_pathway_KEGG.tsv",quote = F, sep = '\t',row.names = F, col.names = T)

enrichment_chart(result_df = bl123_in1235_pathway_KEGG, top_terms = 78)

write.table(bl123_in1235_pathway_KEGG,"/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_in_bl_Enrichment_KEGG.tsv",quote = F, sep = '\t',row.names = F, col.names = T)
bl123_in1235_pathway_GO <- run_pathfindR(
  input = bl123_in1235.pathway,
  convert2alias = FALSE,
  gene_sets = "GO-All",
  # custom_genes = kegg_genes,
  # custom_descriptions = kegg_descriptions,
  iterations = 25,
  n_processes = 8,
  output_dir = "/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_in_bl_GO"
) # Found 372/1036 genes, found 91 enriched GO-All terms
enrichment_chart(result_df = bl123_in1235_pathway_GO, top_terms = 97)

write.table(bl123_in1235_pathway_GO,"/Users/sjblair/Projects/sysAct_allFiles/Mar26_2024_blastemaPairwise/pathfindR_ATAC_in_bl_pathway_GO.tsv",quote = F, sep = '\t',row.names = F, col.names = T)



# I'll make heatmaps now
bl123_cl345.dds@assays$counts
tmp <- bl123_cl345.cts
heatmap_names <- lapply(
  rownames(tmp),
  function(x) bquote(italic(.(x))))
heat <- scales::rescale(t(scale(t(tmp))), to = c(-1,1))
#draw heatmap allowing larger margins and adjusting row label font size
library(pheatmap)
heatmap <- pheatmap(heat,
                    scale="none",
                    cluster_rows=T,
                    cluster_cols=T,
                    show_rownames = T,
                    main = "WAT KI vs WT",
                    color = viridis(25, direction = -1),
                    border_color = NA,
                    labels_row = as.expression(heatmap_names),
                    treeheight_row = 25,
                    treeheight_col =  15)
dev.off()


# I will create a PCA of all of the used replicates
allReplicates.cts <- cts
allReplicates.cts$in4 <- NULL

allReplicates.pca <- prcomp(t(allReplicates.cts), center = T, scale. = T)

# Select the first two principal components
pc1 <- allReplicates.pca$x[, 1]
pc2 <- allReplicates.pca$x[, 2]

# Create the plot
allReplicates.pca <- ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(allReplicates.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC all counts, with Brian Haas") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label=colnames(allReplicates.cts)), size = 5,nudge_x = 5,nudge_y = 10) +
  theme_classic()


# create the new counts from MACS2
setwd("~/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS")

# counts <- read.table("/Users/sjblair/Projects/sysAct_allFiles/ATAC_Mar5_2024/ATAC_COUNTS_Mar15_2024/bl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("Chr","Start","End","Strand","Length","bl1"))
# counts[1:5] <- NULL
bl1 <- read.table("bl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl1"))
bl2 <- read.table("bl2.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl2"))
bl3 <- read.table("bl3.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl3"))
bl4 <- read.table("bl4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl4"))
bl5 <- read.table("bl5.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","bl5"))
cl1 <- read.table("cl1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","cl1"))
cl2 <- read.table("cl2.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","cl2"))
cl3 <- read.table("cl3.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","cl3"))
cl4 <- read.table("cl4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","cl4"))
cl5 <- read.table("cl5.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","cl5"))
in1 <- read.table("in1.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","in1"))
in2 <- read.table("in2.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","in2"))
in3 <- read.table("in3.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","in3"))
in4 <- read.table("in4.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","in4"))
in5 <- read.table("in5.atac.counts", header = T, sep = "\t", row.names = 1, col.names = c("peak","Chr","Start","End","Strand","Length","in5"))

MACS2.cts <- bl1
MACS2.cts$bl2 <- bl2$bl2
MACS2.cts$bl3 <- bl3$bl3
MACS2.cts$bl4 <- bl4$bl4
MACS2.cts$bl5 <- bl5$bl5
MACS2.cts$cl1 <- cl1$cl1
MACS2.cts$cl2 <- cl2$cl2
MACS2.cts$cl3 <- cl3$cl3
MACS2.cts$cl4 <- cl4$cl4
MACS2.cts$cl5 <- cl5$cl5
MACS2.cts$in1 <- in1$in1
MACS2.cts$in2 <- in2$in2
MACS2.cts$in3 <- in3$in3
MACS2.cts$in4 <- in4$in4
MACS2.cts$in5 <- in5$in5
rm(bl1,bl2,bl3,bl4,bl5,in1,in2,in3,in4,in5,cl1,cl2,cl3,cl4,cl5)


pca <- prcomp(t(MACS2.cts[6:20]), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]

# Create the plot
ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(MACS2.cts[6:20])), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC per Wei, 2021") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label =colnames(MACS2.cts[6:20])), size = 5,nudge_x = 2,nudge_y = 2) +
  theme_classic()

# omg, that's ugly, lets remove bl4
MACS2_nobl4.cts <- MACS2.cts
MACS2_nobl4.cts$bl4 <- NULL

pca <- prcomp(t(MACS2_nobl4.cts[6:19]), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]

# Create the plot
ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(MACS2_nobl4.cts[6:19])), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC per Wei, 2021, bl4 omit") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label =colnames(MACS2_nobl4.cts[6:19])), size = 5,nudge_x = 2,nudge_y = 2) +
  theme_classic()

# let's remove the other bad replicates
MACS2_removePoor.rep.cts <- MACS2_nobl4.cts
MACS2_removePoor.rep.cts$cl1 <- NULL
MACS2_removePoor.rep.cts$cl2 <- NULL
MACS2_removePoor.rep.cts$in4 <- NULL
MACS2_removePoor.rep.cts$bl5 <- NULL

MACS2_removePoor.rep.cts$in5 <- NULL
pca <- prcomp(t(MACS2_removePoor.rep.cts), center = TRUE, scale. = T)

# Select the first two principal components
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]

# Create the plot
ggplot(data.frame(PC1 = pc1, PC2 = pc2), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = colnames(MACS2_removePoor.rep.cts)), size = 3) +  # Color points by sample name
  ggtitle("PCA Plot of ATAC per Wei, 2021, bl45,in45,cl12 omit") +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  geom_text_repel(aes(label =colnames(MACS2_removePoor.rep.cts)), size = 5,nudge_x = 2,nudge_y = 2) +
  theme_classic()

# lets try to annotate these peaks
MACS2.anno <- MACS2_removePoor.rep.cts[1:3]
MACS2.anno$peak <- rownames(MACS2.anno)
MACS2.annod <- makeGRangesFromDataFrame(MACS2.anno)

MACS2.anno <- annotatePeak(MACS2.annod,tssRegion=c(-3000, 3000),TxDb=txdb) 
MACS2.anno <- as.data.frame(MACS2.anno)
MACS2.anno$grange <- paste0(MACS2.anno$seqnames,":",MACS2.anno$start,"-",MACS2.anno$end)

MACS2.cbinder <- MACS2_removePoor.rep.cts[1:3]
MACS2.cbinder$peak <- rownames(MACS2.cbinder)
MACS2.cbinder$grange <- paste0(MACS2.cbinder$Chr,":",MACS2.cbinder$Start,"-",MACS2.cbinder$End)


MACS2.bindanno <- left_join(MACS2.anno, MACS2.cbinder)
rownames <- MACS2.bindanno
rownames[1:12] <- NULL
rownames[2:6] <- NULL

rownames$transcriptId <- gsub(" \\[","..",rownames$transcriptId)
rownames$transcriptId <- gsub("\\]","..",rownames$transcriptId)
rownames$transcriptId <- gsub("..nr..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("..hs..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("..hs..","^",rownames$transcriptId)
rownames$transcriptId <- gsub("\\|","^",rownames$transcriptId)
rownames$transcriptId <- gsub("\\^\\^","^",rownames$transcriptId)
rownames$rownames <- paste0(rownames$transcriptId,"^",rownames$peak)

rownames$transcriptId <- NULL
MACS2_removePoor.rep.cts$peak <- rownames(MACS2_removePoor.rep.cts)
MACS2_removePoor.rep.cts <- left_join(MACS2_removePoor.rep.cts,rownames)
MACS2_removePoor.rep.cts$peak <- NULL
rownames(MACS2_removePoor.rep.cts) <- make.names(MACS2_removePoor.rep.cts$rownames, unique = TRUE)
MACS2_removePoor.rep.cts[1:5] <- NULL
MACS2_removePoor.rep.cts$rownames <- NULL
# DA analysis
MACS2_withomisions.cts$in5 <- NULL
# let's see what happens when we remove intact 5
MACS2_withomisions.cts <- MACS2_removePoor.rep.cts
coldata <- colnames(MACS2_withomisions.cts)
condition <- c("Regenerating","Regenerating","Regenerating","Activated","Activated","Activated","Homeostatic","Homeostatic","Homeostatic")
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$coldata
coldata$condition <- condition
names(coldata)[1] <- "sample"
coldata <- as.matrix(coldata)
MACS2_withomisions.dds <- DESeqDataSetFromMatrix(MACS2_withomisions.cts, colData = coldata,
                                                 design = ~ condition)
MACS2_withomisions.dds <- DESeq(MACS2_withomisions.dds)
MACS2_withomisions.sig.act.res <- results(MACS2_withomisions.dds, contrast=c("condition","Regenerating","Activated"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

MACS2_withomisions.sig.act.res <- MACS2_withomisions.sig.act.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)


MACS2_withomisions.sig.regen.res <- results(MACS2_withomisions.dds, contrast=c("condition","Regenerating","Homeostatic"))
# padj_cutoff <- 0.05
# lfc_cutoff <- 0.58 #log(1.5) (base2)

MACS2_withomisions.sig..regenres <- MACS2_withomisions.sig.regen.res %>%
  data.frame()  %>% 
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05)
getwd()

save.image("/Users/sjblair/Projects/sysAct_allFiles/MACS_IDR/MACS2_COUNTS/ATAC.globalEnvironment.Apr09.2024.rdata")




ATF <- c(
  "ATF1-AMEX60DD029673", # ATF-1 
  "ATF1-AMEX60DD035554", # ATF-1 
  "ATF2-AMEX60DD055620", # ATF-2
  "ATF3-AMEX60DDU001003554", # ATF-3
  "ATF4-AMEX60DD029230", # ATF-4
  "ATF5-AMEX60DD017394", # ATF-5
  "ATF6-AMEX60DD018411", # ATF-6
  "ATF6B-AMEX60DD010499", # ATF-6B
  "ATF7-AMEX60DD029531", # ATF-7
  "ATF7IP-AMEX60DD029096", # ATF7IP
  "ATF7IP2-AMEX60DD021479", # ATF7IP2
  "HEATR5B-AMEX60DD011156", # Possible BATF
  "HEATR5B-AMEX60DD011157", # Possible BATF
  "TTLL5-AMEX60DD011157", # Possible BATF
  "BATF2-AMEX60DD025272", # BATF-2
  "BATF3-AMEX60DD036067", # BATF-3
  "JDP2-AMEX60DD011155" # JDP2
)
DotPlot(seu,assay = "RNA",ATF,cols="RdBu",split.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))  + coord_flip()




GENES <- c(
  "EN1-AMEX60DD056025", # En1
  "LHX9-AMEX60DD018728", # Lhx9
  "PITX1-AMEX60DD028957", # Pitx1
  "SMAD3-AMEX60DD003401", # Smad3
  "TWIST2-AMEX60DD055512", # Twist2
  "TWIST3-TWIST2-AMEX60DD029436", # Twist2
  "C-FOS-FOS-AMEX60DD011151", #c-FOS/FOS (AP1)
  "C-FOS-FOS-AMEX60DD011154", #c-FOS/FOS (AP1)
  #"C-FOS-FOS-AMEX60DD011155", #c-FOS
  "ATF3-AMEX60DDU001003554", # ATF-3
  "DLX2-AMEX60DD055671", # DLX2
  "DLX5-AMEX60DD022355", # DLX5
  # "C-FOS-FOS-AMEX60DD011151", #c-FOS
  # "C-FOS-FOS-AMEX60DD011154", #c-FOS/FOS (AP1)
  "FOSL2-AMEX60DD036324", # FOSL2 FRA-2
  "HOXA11-HOXA3-AMEX60DD021947", # Hoxa10
  "JUN-AMEX60DD019448", # Jun-AP1
  "JUN-AMEX60DD043789", # Jun-AP1
  "JUNB-AMEX60DD031925", # JunB
  "NDRG1-AMEX60DD040549", # n-Myc
  "NDRG1-AMEX60DD044627", # n-Myc
  "POU3F1-AMEX60DD005829", #Oct:Oct-short
  "POU2F2-AMEX60DD018055", # Oct2 POU2F2
  "AB205-0109740-PAX6-AMEX60DD007993", # Pax6
  "PAX6-AMEX60DD004905", # Pax6
  "ZNF683-PRDM1-AMEX60DD005654", # Prdm1
  "PRDM1-AMEX60DD034088", # Prdm1
  "SIX2-AMEX60DD035940", # Six2
  "TP53-AMEX60DD012903", # p53 (TP53)
  "WT1-AMEX60DDU001006343", # WT1
  "CLOCK-AMEX60DD045527", # CLOCK
  "CLOCK-AMEX60DD045528", # CLOCK
  "CLOCK-AMEX60DD045529", # CLOCK
  "PASD1-CLOCK-AMEX60DD037299", # CLOCK
  "LIM-1-LHX1-AMEX60DD054124", # LHX1
  "LHX3-AMEX60DD050656", # LHX3
  "POU5F1-AMEX60DD010762", #Oct4-Sox2-Tcf-Nanog
  "POU2-POU5F1-AMEX60DD050901", #Oct4-Sox2-Tcf-Nanog
  "POU2F3-AMEX60DD053808", # OCT11
  #"POU3F1-AMEX60DD005829", # OCT6
  "PAX5-AMEX60DD043936", #PAX5
  "SMAD2-AMEX60DD041720", # SMAD2
  "SMAD4-AMEX60DD015011", # SMAD4
  "SMAD4-AMEX60DD041763", # SMAD4
  #"JUN-AMEX60DD019448", # c-Jun
  #"JUN-AMEX60DD043789", # c-Jun
  "C-MYC-MYC-AMEX60DD040490", # c-Myc
  #"HOXA11-HOXA3-AMEX60DD021947",
  "HOXA13-AMEX60DD021953",
  "HOXB13-AMEX60DD010173", # HOXB13
  # Hoxd12 Not present
  "HOXD13-AMEX60DD055612",
  "LEP-AMEX60DD007977", # Leptin
  "MEIS1-AMEX60DD035598", # Meis1
  "SIX1-AMEX60DD011856", # Six1
  "SIX1-AMEX60DD024753" # Six1
)
GENES <- rev(GENES)
id <- c("EN1","LHX9","PITX1","SMAD3","TWIST2|AMEX60DD055512","TWIST2|AMEX60DD029436","FOS|AMEX60DD011151", "FOS|AMEX60DD011154", "ATF3","DLX2","DLX5","FOSL2","HOXA3","JUN|AMEX60DD019448","JUN|AMEX60DD043789","JUNB","NDRG1|AMEX60DD040549","NDRG1|AMEX60DD044627","POU3F1","POU2F2","PAX6|AMEX60DD007993","PAX6|AMEX60DD004905","PRDM1|AMEX60DD005654","PRDM1|AMEX60DD034088","SIX2","TP53","WT1","CLOCK|AMEX60DD045527","CLOCK|AMEX60DD045528","CLOCK|AMEX60DD045529","CLOCK|AMEX60DD037299","LHX1","LHX3","POU5F1|AMEX60DD010762","POU5F1|AMEX60DD050901","POU2F3|AMEX60DD053808","PAX5|AMEX60DD043936","SMAD2|AMEX60DD041720","SMAD4|AMEX60DD015011","SMAD4|AMEX60DD041763","MYC","HOXA13","HOXB13","HOXD13","LEP","MEIS1","SIX1|AMEX60DD011856","SIX1|AMEX60DD024753")
id <- rev(id)
# Identify cells belonging to the desired clusters (replace "cluster1" and "cluster2" with your actual cluster names)
clusters <- c("Fibroblast I", "Fibroblast II", "Fibroblast III", "Fibroblast IV","Fibroblast V","Fibroblast VI","Fibroblast VII","Fibroblast VIII")

# Create the DotPlot with cluster filtering
DotPlot(seu, assay = "RNA", features = GENES, 
        cols = "RdBu", split.by = "sample", # Remove split.by for single plot
        idents = clusters) +  # Filter cells based on cluster membership
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  coord_flip() +
  scale_x_discrete(labels = id)

clusters <-  c("Epidermal I", "Epidermal II","Epidermal III","Epidermal IV", "Basal epidermis")
DotPlot(seu, assay = "RNA", features = GENES, 
        cols = "RdBu", split.by = "sample", # Remove split.by for single plot
        idents = clusters) +  # Filter cells based on cluster membership
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  coord_flip() +
  scale_x_discrete(labels = id)

levels(seu)

levels(seu) <- c("Fibroblast I","Fibroblast II","Fibroblast III","Fibroblast IV","Fibroblast V","Fibroblast VI","Fibroblast VII","Fibroblast VIII","Epidermal I","Epidermal II","Epidermal III","Epidermal IV","Endothelial I","Endothelial II","Hematopoetic progenitor I","Hematopoetic progenitor II","RBC I","RBC II","Myeloid progenitor","Skeletal muscle","Basal epidermis","T-cell")


clusters <-  c("Hematopoetic progenitor I", "Hematopoetic progenitor II","Endothelial I","Endothelial II", "RBC I", "RBC II", "Myeloid progenitor", "Skeletal muscle", "T-cell")
cluster_order <- factor(clusters, levels = clusters)
DotPlot(seu, assay = "RNA", features = GENES, 
        cols = "RdBu", split.by = "sample", # Remove split.by for single plot
        idents = cluster_order) +  # Filter cells based on cluster membership
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  coord_flip() +
  scale_x_discrete(labels = id)
