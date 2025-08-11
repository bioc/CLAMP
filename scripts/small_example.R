library(PLIER2)

library(rsvd)
library(glmnet)
data("dataWholeBlood")
data("majorCellTypes")

dataWholeBlood=tscale(dataWholeBlood)


gmtList=list()

gmtList[["CellMarkers"]]=getGMT("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024")

gmtList[["KEGG"]]=getGMT("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human", "KEGG_2021_Human")



pathMatCell = gmtListToSparseMat(gmtList)


matchedPathsWB=getMatchedPathwayMat(pathMatCell, new.genes = rownames(dataWholeBlood))

set.seed(1);dataWholeBlood.svd=rsvd(dataWholeBlood,k=25)

matchedPathsWB=readRDS("~/mathedPathsWB.RDS")


wb.plier.old=PLIER::PLIER(dataWholeBlood, as.matrix(matchedPathsWB), svdres = dataWholeBlood.svd, Chat = as.matrix(ChatWB) , doCrossval = T, k=25, pathwaySelection = "fast", maxPath = 10, max.iter=200)




ChatWB=getChat(scale(matchedPathsWB))
wb.plier.base <- PLIERbase(dataWholeBlood, k = 25, svdres = dataWholeBlood.svd, trace = T, adaptive.p = 0.05)

wb.plier.full=PLIERfull(dataWholeBlood,priorMat =  matchedPathsWB, plier.base.result  = wb.plier.base,
                        Chat = ChatWB,doCrossval = T, max.U.updates = 50, trace=T)
celltypeTargets=as.matrix((as.data.frame(majorCellTypes)[colnames(dataWholeBlood),]))

#output=compareBs(wb.plier.base, wb.plier.full,  celltypeTargets, "PLIERv1", "PLIERv2", m="p")
#output$plot

library(ggplot2)
library(ggrepel)
output=compareBs(wb.plier.old, wb.plier.full,  celltypeTargets, "PLIERv1", "PLIERv2", m="p")
output$plot

plotTopZ_Complex(wb.plier.full, data = dataWholeBlood, priorMat = matchedPathsWB, index = c(16,17,24), top = 10, Zheat = T)
