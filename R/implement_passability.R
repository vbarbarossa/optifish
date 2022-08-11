# model implementation
library(igraph); library(dplyr)

HYBAS <- data.frame(INTER_ID = c(1,2,3,4,5), INTER_NEXT = c(0,1,2,1,4))
l = c(20,10,5,5,10)
p = c(1,0.5,0.3,0.6,0.3)

# Cote et al example
HYBAS <- data.frame(INTER_ID = c(1,2), INTER_NEXT = c(0,1))
l = c(20,10)
p = c(1,0.5)

df <- data.frame(from = HYBAS$INTER_ID,to = HYBAS$INTER_NEXT)
outlet <- as.character(df$from[which(df$to == 0)])

df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g <- graph_from_data_frame(d=df)
plot(g)
E(g)$pass <- p

vertices_id <- names(igraph::V(g))
C <- 10^dodgr::dodgr_dists(
  mutate(
    rbind(
      igraph::as_data_frame(igraph::as.undirected(g,mode = 'each'), what = "edges")
      ,
      select(
        rename(
          igraph::as_data_frame(igraph::as.undirected(g,mode = 'each'), what = "edges"),
          from = to, to = from
        ),from,to,pass
      )
    ),dist = log10(.data$pass)
  ),
  from = vertices_id, to = vertices_id)

I = matrix(1, nrow = length(l), ncol = length(l))
L = (I*l * t(I*l)) / sum(l)**2 

sum(C*L)*100

# diadromous
sum(C[,1]*(l/sum(l)))*100

# avoid using dodgr
# pw <- -1*ifelse(p!=0, log(p), 0)
pw <- -1*log(p)
C2 <- exp(distances(g, mode = 'all', weights = pw) * -1)
round(c(C),3) == round(c(C2),3)
# yes!

