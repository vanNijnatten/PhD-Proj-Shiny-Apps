library(shiny)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(pheatmap)
library(Seurat)
library(DESeq2)
library(RColorBrewer)
library(readxl)
library(ggsignif)
library(DT)
library(graphics)
library(enrichplot)
library(ggpubr)
library(ggprism)
library(patchwork)

##Normalised Counts file
norm<-read.csv("data/normalizedcounts_CPM.csv")
rownames(norm)<-norm[,1]
norm<-norm[,-1]
exp<-norm

##DEG results - for Volcano plots and DEG table
ava<-as.data.frame(read_excel("data/fullDEG_allcomparisons_newFC.xlsx",
                sheet="airgsk547_vs_airvehicle"))
rownames(ava)<-ava[,1]
ava<-ava[,-1]

avcs_gsk<-as.data.frame(read_excel("data/fullDEG_allcomparisons_newFC.xlsx",
                     sheet="airgsk547_vs_csgsk547"))
rownames(avcs_gsk)<-avcs_gsk[,1]
avcs_gsk<-avcs_gsk[,-1]

avcs_veh<-as.data.frame(read_excel("data/fullDEG_allcomparisons_newFC.xlsx",
                     sheet="airvehicle_vs_csvehicle"))
rownames(avcs_veh)<-avcs_veh[,1]
avcs_veh<-avcs_veh[,-1]

cvc<-as.data.frame(read_excel("data/fullDEG_allcomparisons_newFC.xlsx",
                sheet="csgsk547_vs_csvehicle"))
rownames(cvc)<-cvc[,1]
cvc<-cvc[,-1]


##GSEA results
go_active<-readRDS("data/GSEA_go_activated.rds")
go_supressed<-readRDS("data/GSEA_go_supressed.rds")

kegg_active<-readRDS("data/GSEA_kegg_activated.rds")
kegg_supressed<-readRDS("data/GSEA_kegg_supressed.rds")





ui<-fluidPage(

  titlePanel("CD29 Project"),

  sidebarLayout(

      sidebarPanel(

        h1("Differential Expression Results"),
          selectInput(inputId = "comparisons",
                      label="Comparisons",
                      choices=c("Air GSK547 vs Air Vehicle",
                                "Air GSK547 vs CS GSK547",
                                "Air Vehicle vs CS Vehicle",
                                "CS GSK547 vs CS Vehicle"),
                      selected = NULL,
                      multiple = FALSE),
          radioButtons(inputId = "label",
                       label="Volcano Plots Label",
                       choices=c("Show","Hide"),
                       selected=NULL),
          radioButtons(inputId="enrich",
                       label="GSEA",
                       choices=c("Activated","Supressed"),
                       select=NULL),
        radioButtons(inputId = "order",
                     label="GSEA ordering",
                     choices=c("P-Val Adjusted","Normalised Enrichment Score"),
                     select=NULL),
        selectInput(inputId ="exp",
                    label="Genes",
                    choices=rownames(exp),
                    selected="Lypla1",
                    multiple = FALSE)
    ),

    mainPanel(tabsetPanel(
      tabPanel("Volcano Plots",
               plotOutput("volcano")),
      tabPanel("Differential Expression Results Table",
               DT::dataTableOutput("deg")),
      tabPanel("GSEA - GO",
               plotOutput("gsea_go")),
      tabPanel("GSEA - KEGG",
               plotOutput("gsea_kegg")),
      tabPanel("Normalised Gene Expression",
               plotOutput("exp"))
      ))

    )

)



server<-function(input,output){

  output$volcano<-renderPlot({

    if(input$comparisons=="Air GSK547 vs Air Vehicle"){
      p1=ggplot(ava, aes(logFC, -log10(FDR))) +
        geom_point(aes(colour = Legend)) +
        scale_color_manual(values = c("#07BDF4","#000000","#FB031C")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed") +
        geom_vline(xintercept = -0.58, linetype="dashed") +
        geom_vline(xintercept = 0.58, linetype="dashed") +
        xlab(expression(bold(Log[2]*"FC")))+
        ylab(expression(bold(-Log[10]*"(FDR)")))+
        theme(panel.background=element_rect(fill="white"),
              panel.grid=element_line(colour="white", size=0.15),
              #plot.title=element_text(size=17, hjust=0.5),
              axis.title=element_text(size=18,colour="black", face="bold"),
              axis.text = element_text(size=20),
              aspect.ratio = 0.5,
              axis.line.x = element_line(colour="black", size = 0.5),
              axis.line.y = element_line(colour="black", size = 0.5))+
        scale_x_continuous(minor_breaks = waiver(), breaks = c(seq(-8, 8, by= 2)))+
        scale_y_continuous(minor_breaks = FALSE)+
        coord_cartesian( ylim = c(0,5),xlim=c(-5,5),clip="off")+NoLegend()

      if(input$label=="Show"){
        resultsf <- ava[which((ava$logFC>0.58|ava$logFC<(-0.58))&
                                    ava$FDR<0.05),]
        tT2<-resultsf
       print(p1+ geom_text_repel(data = tT2[head(order(tT2$logFC), n=5),],
          aes(label=rownames(tT2[head(order(tT2$logFC), n=5),])),
          min.segment.length = 0, nudge_x = -0.4, direction = "y", hjust = "right")+
          geom_text_repel(data = tT2[tail(order(tT2$logFC), n=5),],
          aes(label=rownames(tT2[tail(order(tT2$logFC), n=5),])),
          min.segment.length = 0,nudge_x = 0.4, direction = "y", hjust = "left"))
      }
      if(input$label=="Hide"){
        print(p1)
      }
    }
    if(input$comparisons=="Air GSK547 vs CS GSK547"){

      p2=ggplot(avcs_gsk, aes(logFC, -log10(FDR))) +
        geom_point(aes(colour = Legend)) +
        scale_color_manual(values = c("#07BDF4","#000000","#FB031C")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed") +
        geom_vline(xintercept = -0.58, linetype="dashed") +
        geom_vline(xintercept = 0.58, linetype="dashed") +
        xlab(expression(bold(Log[2]*"FC")))+
        ylab(expression(bold(-Log[10]*"(FDR)")))+
        theme(panel.background=element_rect(fill="white"),
              panel.grid=element_line(colour="white", size=0.15),
              #plot.title=element_text(size=17, hjust=0.5),
              axis.title=element_text(size=18,colour="black", face="bold"),
              axis.text = element_text(size=20),
              aspect.ratio = 0.5,
              axis.line.x = element_line(colour="black", size = 0.5),
              axis.line.y = element_line(colour="black", size = 0.5))+
        scale_x_continuous(minor_breaks = waiver(), breaks = c(seq(-8, 8, by= 2)))+
        scale_y_continuous(minor_breaks = FALSE)+
        coord_cartesian( ylim = c(0,20),xlim=c(-9,9),clip="off")+NoLegend()

      if(input$label=="Show"){
        resultsf <- avcs_gsk[which((avcs_gsk$logFC>0.58|avcs_gsk$logFC<(-0.58))&
                                avcs_gsk$FDR<0.05),]
        tT2<-resultsf
        print(p2+ geom_text_repel(data = tT2[head(order(tT2$logFC), n=5),],
                                  aes(label=rownames(tT2[head(order(tT2$logFC), n=5),])),
                                  min.segment.length = 0, nudge_x = -0.4, direction = "y", hjust = "right")+
                geom_text_repel(data = tT2[tail(order(tT2$logFC), n=5),],
                                aes(label=rownames(tT2[tail(order(tT2$logFC), n=5),])),
                                min.segment.length = 0,nudge_x = 0.4, direction = "y", hjust = "left"))
      }
      if(input$label=="Hide"){
        print(p2)
      }
    }

    if(input$comparisons=="Air Vehicle vs CS Vehicle"){
     p3= ggplot(avcs_veh, aes(logFC, -log10(FDR))) +
        geom_point(aes(colour = Legend)) +
        scale_color_manual(values = c("#07BDF4","#000000","#FB031C")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed") +
        geom_vline(xintercept = -0.58, linetype="dashed") +
        geom_vline(xintercept = 0.58, linetype="dashed") +
        xlab(expression(bold(Log[2]*"FC")))+
        ylab(expression(bold(-Log[10]*"(FDR)")))+
        theme(panel.background=element_rect(fill="white"),
              panel.grid=element_line(colour="white", size=0.15),
              #plot.title=element_text(size=17, hjust=0.5),
              axis.title=element_text(size=18,colour="black", face="bold"),
              axis.text = element_text(size=20),
              aspect.ratio = 0.5,
              axis.line.x = element_line(colour="black", size = 0.5),
              axis.line.y = element_line(colour="black", size = 0.5))+
        scale_x_continuous(minor_breaks = waiver(), breaks = c(seq(-8, 8, by= 2)))+
        scale_y_continuous(minor_breaks = FALSE)+
        coord_cartesian( ylim = c(0,25),xlim=c(-8,8),clip="off")+NoLegend()

     if(input$label=="Show"){
       resultsf <- avcs_veh[which((avcs_veh$logFC>0.58|avcs_veh$logFC<(-0.58))&
                                    avcs_veh$FDR<0.05),]
       tT2<-resultsf
       print(p3+ geom_text_repel(data = tT2[head(order(tT2$logFC), n=5),],
                                 aes(label=rownames(tT2[head(order(tT2$logFC), n=5),])),
                                 min.segment.length = 0, nudge_x = -0.4, direction = "y", hjust = "right")+
               geom_text_repel(data = tT2[tail(order(tT2$logFC), n=5),],
                               aes(label=rownames(tT2[tail(order(tT2$logFC), n=5),])),
                               min.segment.length = 0,nudge_x = 0.4, direction = "y", hjust = "left"))
     }
     if(input$label=="Hide"){
       print(p3)
     }
  }

    if(input$comparisons=="CS GSK547 vs CS Vehicle"){
      p4=ggplot(cvc, aes(logFC, -log10(FDR))) +
        geom_point(aes(colour = Legend)) +
        scale_color_manual(values = c("#07BDF4","#000000","#FB031C")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed") +
        geom_vline(xintercept = -0.58, linetype="dashed") +
        geom_vline(xintercept = 0.58, linetype="dashed") +
        xlab(expression(bold(Log[2]*"FC")))+
        ylab(expression(bold(-Log[10]*"(FDR)")))+
        theme(panel.background=element_rect(fill="white"),
              panel.grid=element_line(colour="white", size=0.15),
              #plot.title=element_text(size=17, hjust=0.5),
              axis.title=element_text(size=18,colour="black", face="bold"),
              axis.text = element_text(size=20),
              aspect.ratio = 0.5,
              axis.line.x = element_line(colour="black", size = 0.5),
              axis.line.y = element_line(colour="black", size = 0.5))+
        scale_x_continuous(minor_breaks = waiver(), breaks = c(seq(-8, 8, by= 2)))+
        scale_y_continuous(minor_breaks = FALSE)+
        coord_cartesian( ylim = c(0,25),xlim=c(-4,4),clip="off")+NoLegend()

      if(input$label=="Show"){
        resultsf <- cvc[which((cvc$logFC>0.58|cvc$logFC<(-0.58))&
                                     cvc$FDR<0.05),]
        tT2<-resultsf
        print(p4+ geom_text_repel(data = tT2[head(order(tT2$logFC), n=5),],
                                  aes(label=rownames(tT2[head(order(tT2$logFC), n=5),])),
                                  min.segment.length = 0, nudge_x = -0.4, direction = "y", hjust = "right")+
                geom_text_repel(data = tT2[tail(order(tT2$logFC), n=5),],
                                aes(label=rownames(tT2[tail(order(tT2$logFC), n=5),])),
                                min.segment.length = 0,nudge_x = 0.4, direction = "y", hjust = "left"))
      }
      if(input$label=="Hide"){
        print(p4)
      }
  }

  })

  output$deg<-DT::renderDataTable({
    if(input$comparisons=="Air GSK547 vs Air Vehicle"){
      datatable(ava,rownames=TRUE)
    }
    else if(input$comparisons=="Air GSK547 vs CS GSK547"){
      datatable(avcs_gsk,rownames=TRUE)
    }
    else if(input$comparisons=="Air Vehicle vs CS Vehicle"){
      datatable(avcs_veh,rownames=TRUE)
    }
    else if(input$comparisons=="CS GSK547 vs CS Vehicle"){
      datatable(cvc,rownames = TRUE)
    }
   })

  output$exp<-renderPlot({

    d1<-data.frame(Group=c(rep("Air GSK547",6),
                           rep("Air Vehicle",6),
                           rep("CS GSK547",6),
                           rep("CS Vehicle",6)),
                   t(norm))

   ggboxplot(d1,x="Group",y=input$exp,
             fill=c("#D1E5F0","salmon","#2166AC","#b2182b"))+
     theme(axis.text.x = element_text(color = "black", size = 16, face = "bold"),
           axis.text.y = element_text(color = "black", size = 16, angle = 0,
                                      hjust = 1, vjust = 0, face = "plain"),
           axis.title.x = element_blank(),
           axis.title.y = element_text(color = "black", size = 16, angle = 90,
                                       hjust = .5, vjust = 2, face = "bold"),
           legend.position = "none")+
     scale_y_continuous(name= input$exp)+ #change
     scale_x_discrete(limits = c("Air GSK547","Air Vehicle",
                                 "CS GSK547","CS Vehicle"))+
     geom_bracket(xmin="Air GSK547",xmax="Air Vehicle",y.position=max(d1[,input$exp])+5,
                  label=signif(ava[input$exp,"PValue"],digits=3))+
     geom_bracket(xmin="Air GSK547",xmax="CS GSK547",y.position=max(d1[,input$exp])+10,
                  label=signif(avcs_gsk[input$exp,"PValue"],digits=3))+
     geom_bracket(xmin="Air Vehicle",xmax="CS Vehicle",y.position=max(d1[,input$exp])+18,
                  label=signif(avcs_veh[input$exp,"PValue"],digits=3))+
     geom_bracket(xmin="CS GSK547",xmax="CS Vehicle",y.position=max(d1[,input$exp])+5,
                  label=signif(cvc[input$exp,"PValue"],digits=3))
  })

  output$gsea_go<-renderPlot({

    if(input$comparisons=="Air GSK547 vs Air Vehicle"){
      if(input$enrich=="Activated"){
      active <- dotplot(go_active[[1]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
        facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

      print(active)

      if(input$order=="Normalised Enrichment Score"){
        active <- dotplot(go_active[[1]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
        facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
        print(active)
      }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(go_supressed[[1]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(go_supressed[[1]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="Air GSK547 vs CS GSK547"){
      if(input$enrich=="Activated"){
        active <- dotplot(go_active[[2]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(go_active[[2]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(go_supressed[[2]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(go_supressed[[2]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="Air Vehicle vs CS Vehicle"){
      if(input$enrich=="Activated"){
        active <- dotplot(go_active[[3]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(go_active[[3]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(go_supressed[[3]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(go_supressed[[3]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="CS GSK547 vs CS Vehicle"){
      if(input$enrich=="Activated"){
        active <- dotplot(go_active[[4]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(go_active[[4]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(go_supressed[[4]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(go_supressed[[4]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }

  })

  output$gsea_kegg<-renderPlot({

    if(input$comparisons=="Air GSK547 vs Air Vehicle"){
      if(input$enrich=="Activated"){
        active <- dotplot(kegg_active[[1]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(kegg_active[[1]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(kegg_supressed[[1]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(kegg_supressed[[1]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="Air GSK547 vs CS GSK547"){
      if(input$enrich=="Activated"){
        active <- dotplot(kegg_active[[2]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(kegg_active[[2]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(kegg_supressed[[2]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(kegg_supressed[[2]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="Air Vehicle vs CS Vehicle"){
      if(input$enrich=="Activated"){
        active <- dotplot(kegg_active[[3]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(kegg_active[[3]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(kegg_supressed[[3]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(kegg_supressed[[3]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }
    if(input$comparisons=="CS GSK547 vs CS Vehicle"){
      if(input$enrich=="Activated"){
        active <- dotplot(kegg_active[[4]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(active)

        if(input$order=="Normalised Enrichment Score"){
          active <- dotplot(kegg_active[[4]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(active)
        }
      }
      if(input$enrich=="Supressed"){
        supressed <- dotplot(kegg_supressed[[4]],showCategory=10,split=".sign",orderBy="Description",x="p.adjust") +
          facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")

        print(supressed)

        if(input$order=="Normalised Enrichment Score"){
          supressed <- dotplot(kegg_supressed[[4]],showCategory=10,split=".sign",orderBy="Description",x="NES") +
            facet_grid(.~.sign,scale="free")+scale_color_gradientn(colours="red")
          print(supressed)
        }
      }
    }

  })
}


shinyApp(ui,server)
