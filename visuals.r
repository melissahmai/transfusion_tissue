library(dplyr) # Data manipulation
library(ggplot2) # General graphing
library(igraph)
library(qgraph)
library(ggraph) # Graph visualization
library(cowplot)
library(grid) # Plot rearrangement
library(magick)
library(pdftools) # Image manipulation
library(latex2exp) # Latex
library(abind) # N-dimensional concatenation
library(extrafont)

##### Rotations
rotate <- function(layout, theta) {
  trial <- layout
  
  trial$x <- layout$x*cos(theta) - layout$y*sin(theta)
  trial$y <- layout$x*sin(theta) + layout$y*cos(theta)
  
  return(trial)
}


##### Create a network layout object from [.]_nodes.csv and [.]_edges.csv files
net_layout <- function(system, seed = NULL, set.coords = NULL) {
  # Read in the node information
  nodes <- read.csv(paste0(system,'_nodes.csv'))
  
  nodes$Label[nodes$Label=='a'] <- 'air'
  nodes$Label[nodes$Label=='mc'] <- 'ms'
  # Translate cell codes to descriptive classes and locations
  nodes <- nodes %>% mutate(Class = ifelse(Label %in% c('bx','ax','tt'), 
                                           'Tracheary', 
                                           ifelse(Label %in% c('bp','ap','tp'),
                                                  'Sugar',
                                                  ifelse(Label %in% c('bs','ms','mc'),
                                                         'Extrastelar',
                                                         'Air'))),
                            Location = ifelse(Label %in% c('bp','bx','a','ms','mc'),
                                              'Boundary',
                                              ifelse(Label %in% c('ax','ap'),
                                                     'Vascular Bundle',
                                                     ifelse(Label %in% c('tp','tt'),
                                                            'Transfusion Tissue',
                                                            'Bundle Sheath'))))
  
  #
  nodes$Psi[nodes$Label %in% c('bx','bp','air')] <- NA
  
  # Read in edge information
  edges <- read.csv(paste0(system,'_edges.csv'))
  
  # Rename columns (Source and Target are needed for Gephi)
  edges <- edges %>% select(Source,Target,value,flow) %>% rename(from=Source, to=Target)
  
  # Normalize the flows according to type
  edges <- edges %>% group_by(flow) %>% mutate(value=value/max(value,na.rm=T)) %>%
    ungroup() %>% as.data.frame
  
  # Assign a descriptive name to the flow type
  edges$flow <- ifelse(edges$flow=='u','Water',ifelse(edges$flow=='j','Sugar','Contrast'))
  
  # Create the network object 
  net <- graph_from_data_frame(d=edges,vertices=nodes,directed=T)
  
  # Create the layout from the object, using a seed if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  layout <- create_layout(net, layout = 'fr')
  
  # If already defined coordinates are given, overwrite
  if (!is.null(set.coords)) {
    coords <- read.csv(set.coords, header=T)
    coords <- coords[match(layout$name,coords$ID),]
    layout[,c('x','y')] <- coords[,c('x','y')]
  }
  
  return(layout)
}

colormap <- c(
  rgb(0,0,1), #ax
  rgb(0.8,0,0.5), #ap
  rgb(0,0.4,0), #ms
  rgb(0,1,1), #tt
  rgb(0.85,0.45,0), #tp
  rgb(0,0.8,0.2), #bs
  rgb(0,0,0.8),
  rgb(0.64,0,0.4),
  rgb(0.7,0.7,0.7)
)

##### Make a preliminary graph
test_layout <- function(layout, show.labels = F, show.edges = T, edge_size = 1,
                        node_size = 4, show.grid = F, show.fig = T, 
                        limx=NULL, limy=NULL, full.labels = F,
                        bg_color = 'transparent', bold.text = F) {
  
  # Remove airspace water potential to avoid dragging down the colorbar
  layout$Psi[layout$Label=='air'] <- max(layout$Psi)
  
  # Initialize the graph using the inputted layout (from net_layout)
  main <- ggraph(layout)
  
  # Print out the edges underneath the nodes, if displaying
  if (show.edges) {
    main <- main + 
      geom_edge_fan(
        aes(color = flow), end_cap = circle(node_size*2.5, 'pt'), strength = 2, edge_width = edge_size,
        arrow = arrow(type='closed', angle = 15, length = unit(4,'pt'))
      )
  }
  
  if (full.labels) {
    labels <- c('Axial Xylem (ax)', 'Axial Phloem (ap)', 'Mesophyll (ms)', 
                'Transfusion Tracheid (tt)', 'Transfusion Parenchyma (tp)', 
                'Bundle Sheath (bs)',
                'Boundary Xylem (bx)', 'Boundary Phloem (bp)', 'Air')
  } else {
    labels <- c('ax','ap','ms','tt','tp','bs','b','c','a')
  }
  
  shapes <- c(24,25,23,23,21,22,24,24,24) # to be filled
  # shapes <- c(19,19,17,18,18,15,17,17,17) # no filling (no psi)
  main <- main + 
    # Aesthetics
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill = 'transparent', color = 'transparent'),
          legend.background = element_rect(fill = 'transparent', color=NA),
          legend.box.background = element_rect(fill = 'transparent', color=NA),
          legend.key=element_blank(), 
          legend.position = 'bottom',
          legend.direction= 'horizontal', 
          legend.box = 'horizontal') +
    # Nodes
    geom_node_point(aes(color=Label,fill=Label,shape=Label,size=(Class=='Air')),stroke=2) +
    # Legends
    scale_size_manual(values = c(node_size,0),guide = 'none') +
    scale_shape_manual(values = shapes,
                       breaks = c('ax','ap','ms','tt','tp','bs','b','c','a'),#c('ax','tt','bx','ap','tp','bp','ms','bs','air'),
                       labels = labels,
                       guide = guide_legend(title = 'Cell/Node Type', title.position = 'top', nrow = 3)) + 
    scale_color_manual(
      values = colormap,
      breaks = c('ax','ap','ms','tt','tp','bs','b','c','a'),#c('ax','tt','b','ap','tp','c','ms','bs','a'),
      labels = labels,
      guide = guide_legend(title = 'Cell/Node Type', title.position = 'top', nrow = 3)) +
    scale_fill_manual(
      values = colormap,
      breaks = c('ax','ap','ms','tt','tp','bs','b','c','a'),#c('ax','tt','b','ap','tp','c','ms','bs','a'),
      labels = labels,
      guide = guide_legend(title = 'Cell/Node Type', title.position = 'top', nrow = 3))+ 
    # Coordinate limits
    coord_fixed(ratio = 1, xlim = limx, ylim=limy, clip = 'on', expand = T)
  
  # Invert the label colors if background is not white or transparent
  if ((bg_color != 'transparent') & (bg_color != 'white')) {
    main <- main + 
      theme(
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white')
      )
  }
  
  # Showing labels
  if (show.labels) {
    main <- main + geom_node_label(aes(label = name))
  }
  
  # Showing the grid (if, for example, needing to adjust coordinates/find zoom limits)
  if (show.grid) {
    main <- main + theme(axis.line = element_line(color = 'black'),
                         axis.ticks = element_line(color = 'black'),
                         axis.text = element_text(color='black'),
                         panel.grid = element_line(color = '#bbbbbb'))
  }
  
  # If bold text
  if (bold.text) {
     main <- main + theme(
        legend.text = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
     )
  }
  
  # Print out the figure if desired
  if (show.fig) {
    print(main)
  }
  return(main)
}

# Colormap for the water flows
# watermap <- c(
#   hsv(0.45,.125,1),
#   hsv(0.45,.25,1),
#   hsv(0.45,.5,1),
#   hsv(0.45,1,1),
#   hsv(0.475,1,1),
#   hsv(0.5,1,1),
#   hsv(0.525,1,1),
#   hsv(0.55,1,1),
#   hsv(0.6,1,1),
#   hsv(0.65,1,1),
#   hsv(.7,1,1)
# )

watermap <- c(
  rgb(0.8,1,1),
  rgb(0,1,1),
  rgb(0,0,1)
)

# Colormap for the sugar flows
# sugarmap <- c(
#   hsv(0.1,1,1),
#   hsv(0.05,1,1),
#   hsv(0,1,1),
#   hsv(0.9,1,1),
#   hsv(0.85,1,1)
# )
sugarmap <- c(
  rgb(1,1,0.8),
  rgb(1,0.5,0),
  rgb(1,0.2,0.5)
)

# Displaying only the edges
# Input 'flow' should be either 'Water' or 'Sugar'
edge_only <- function(layout, Flow, edge_size = 1, arrow_size = 1, limx=NULL, limy=NULL,
                      node_size = 4, bg_color = 'transparent', bold.text = F) {
  # Initialize the graph with the inputted layout (from net_layout)
  p <- ggraph(layout)
  
  # Choose the right flow
  if (Flow == 'Water') {
    colormap <- watermap
    if (bg_color == 'black') {
      # Reverse the colors for better readability if the background is black
      colormap <- colormap[length(colormap):1]
    }
    # p <- p + geom_edge_fan(aes(alpha = (flow=='Water')),#, color=log10(value)), 
    #                        end_cap = circle(10,'pt'), strength = 1.5, edge_width = edge_size,
    #                        arrow = arrow(type='closed', angle = 30, length = unit(arrow_size*4,'pt')))
  } else {
    colormap = sugarmap
  }
  p <- p + 
    geom_node_point(aes(size=(Class=='Air')), alpha = 0, stroke=2, show.legend = F) +
    geom_edge_fan(aes(alpha = (flow==Flow), color = log10(value)), 
                         end_cap = circle(node_size*2.5,'pt'), strength = 1.5, edge_width = edge_size,
                         arrow = arrow(type='closed', angle = 30, length = unit(arrow_size*4,'pt')))
  
  # Legends and other aesthetics
  p <- p +
    scale_edge_alpha_discrete(range = c(0,1), guide = 'none') + 
    scale_edge_color_gradientn(
      colors = colormap,
      limits = c(-3.5,0),
      guide = guide_edge_colorbar(title = TeX(ifelse(Flow=='Water', 
                                                     '$\\log(U/U_x)$',
                                                     '$\\log(J/J_p)$')),
                                  direction = 'vertical')
    ) + 
    theme(panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = 'transparent'),
          legend.background = element_rect(fill = 'transparent', color=NA),
          legend.box.background = element_rect(fill = 'transparent', color=NA),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          legend.box = 'horizontal',
          legend.justification = c(1,0.5))+ 
    coord_fixed(ratio = 1, xlim = limx, ylim=limy, clip = 'off', expand = T) 
  
  # If bold text
  if (bold.text) {
     p <- p + theme(
        legend.text = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
     )
  }
  
  # Invert label color if background is dark
  if ((bg_color != 'transparent') & (bg_color != 'white')) {
    p <- p + 
      theme(
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white')
      ) 
  }
  
  return(p)
}


# To create the actual figure
makefig <- function(layout, limx=NULL, limy=NULL, 
                    edge_size = 1, arrow_size = 1, node_size = 4,
                    width = NA, height = NA, bg_color = 'transparent',
                    is.transparent = T, bold.text = F, full.labels = F,
                    ext = 'pdf', include.legends = T) {
  
  # Get nodes
  main <- test_layout(layout, show.edges = F, show.fig = F, limx=limx, node_size = node_size,
                      limy=limy, bg_color = bg_color, bold.text=bold.text, full.labels = full.labels) + 
    theme(text=element_text(family='serif'))
  
  # Extract the legend for node features
  mainleg <- cowplot::get_legend(main)
  
  # Remove legends from the  main section
  main <- main + guides(color = 'none', fill = 'none', shape = 'none')
  
  # Get water flows
  water <- edge_only(layout, 'Water', edge_size = edge_size, node_size = node_size,
                     arrow_size = arrow_size, limx=limx, limy=limy,
                     bg_color = bg_color, bold.text=bold.text)  + 
    theme(text=element_text(family='serif'))
  
  # Extract the water flow legend and remove from figure
  wbar <- cowplot::get_legend(water)
  water <- water + guides(edge_color = 'none')
  
  # Get sugar flows
  sugar <- edge_only(layout, 'Sugar', edge_size = edge_size, node_size = node_size,
                     arrow_size = arrow_size, limx=limx, limy=limy,
                     bg_color = bg_color, bold.text=bold.text)  + 
    theme(text=element_text(family='serif'))
  
  # Extract the sguar flow legend and remove from figure
  sbar <- cowplot::get_legend(sugar)
  sugar <- sugar + guides(edge_color = 'none')
  
  # Temporarily save the figures (because magick doesn't do well with converting directly from ggplot?)
  ggsave(paste0('temp/main.',ext), main, width = width, height = height, units = 'px', bg = bg_color,dpi=600)
  ggsave(paste0('temp/sugar.',ext), sugar, width = width, height = height, units = 'px', bg = bg_color,dpi=600)
  ggsave(paste0('temp/water.',ext), water,width = width, height = height, units = 'px', bg = bg_color,dpi=600)
  
  if (include.legends) {
    ggsave(paste0('temp/sbar.',ext), sbar, bg = bg_color, width = width, height = height, units = 'px',dpi=600)
    ggsave(paste0('temp/wbar.',ext), wbar, bg = bg_color, width = width, height = height, units = 'px',dpi=600)
    ggsave('temp/mainleg.pdf', mainleg, bg = bg_color, width = width, height = height, units = 'px',dpi=600)
  }
  
  # Re-read in plots as images and set the transparency
  bg_color = ifelse(bg_color=='transparent','white',bg_color)
  readfun <- function(file,ext) {
     if (ext=='pdf') {
        im <- image_read_pdf(paste0(file,'.',ext))
     }
     if (ext=='svg') {
        im <- image_read_svg(paste0(file,'.',ext))
     }
     
     return(im)
  }
  main <- readfun('temp/main',ext) %>% image_transparent(color = bg_color)
  sugar <- readfun('temp/sugar',ext) %>% image_transparent(color = bg_color)
  water <- readfun('temp/water',ext) %>% image_transparent(color = bg_color)
  
  # Overlay the main part of the figure (full network)
  network <- image_mosaic(c(water,sugar,main)) %>% image_transparent(color = 'white')
  
  if (include.legends) {
    
    # Sidebar legends (sugar and water)
    # Read in the sugar colorbar
    sbar <- readfun('temp/sbar',ext) %>% image_transparent(color = bg_color)
    
    # Get an initial object for a blank filler
    filler <- image_crop(sbar,'x1')
    
    # Remove whitespace around the colorbar
    sbar <- image_trim(sbar)
    
    # Get dimensions of the the trimmed sugar colorbar
    sw <- dim(sbar[[1]])[2]
    
    # Read in the water color bar and trim
    wbar <- readfun('temp/wbar',ext) %>% image_transparent(color = bg_color) %>% image_trim
    
    # Get water colorbar dimensions
    ww <- dim(wbar[[1]])[2]
    
    # The total width of the sidebar is that of the widest colorbar
    barw <- max(sw,ww)
    
    # Crop the filler to be the width of the sidebar
    filler <- image_crop(filler,paste0(barw,'x1'))
    
    # Add filler to the side of the smaller legend
    if (sw > ww) {
      wbar <- image_append(
        c(wbar,image_extent(
          image_crop(filler,'1x1'), 
          paste0(sw-ww,'x',dim(wbar[[1]])[3])))
      )
    } else {
      sbar <- image_append(
        c(sbar,image_extent(
          image_crop(filler,'1x1'), 
          paste0(ww-sw,'x',dim(sbar[[1]])[3])))
      )
    }
    
    # Stack the legends together with some space in between
    sidelegends <- image_append(c(wbar,rep(filler,100),sbar), stack = TRUE) %>% image_trim
    
    # Get full dimensions of the sidebar
    sidedims <- dim(sidelegends[[1]])
    sw <- sidedims[2]
    sh <- sidedims[3]
    
    # Now create a filler the width of the sidebar to center it vertically with the main figure
    filler <- image_crop(filler, paste0(sw,'x'))
    to.fill <- max(0,floor((dim(network[[1]])[3]-sh)/2))
    end.fill <- max(0,dim(network[[1]])[3] - sh - to.fill)
    
    # Add filler above and below to center it with the main figure
    sidelegends <- image_append(c(rep(filler,to.fill),sidelegends,rep(filler,end.fill)),stack=TRUE)
    
    # Append the sidebar to the main figure and get new dimensions
    network <- image_append(c(network,sidelegends))
    netsidedims <- dim(network[[1]])
    nw <- netsidedims[2]
    
    # Read in the node legend and create the filler, as before
    mainleg <- image_read_pdf('temp/mainleg.pdf') %>% image_transparent(color = bg_color)
    filler <- image_crop(mainleg,'1x')
    mainleg <- mainleg %>% image_trim
    filler <- image_crop(filler,paste0('x',dim(mainleg[[1]])[3]))
    to.fill <- floor((nw-dim(mainleg[[1]])[2])/2)
    end.fill <- nw - to.fill - dim(mainleg[[1]])[2]
    
    # Stack the (main+sidebar) with the node legend
    if (to.fill<0) {
      to.fill <- -to.fill
      end.fill <- -end.fill
      network <- image_append(c(image_append(c(rep(filler,to.fill),network,rep(filler,end.fill))),mainleg), stack = TRUE)
    } else {
      network <- image_append(c(network,image_append(c(rep(filler,to.fill),mainleg,rep(filler,end.fill)))), stack = TRUE)
    }
  }
  
  if (!is.transparent) {
    network <- image_background(network,bg_color)
  }
  
  return(network)
}

