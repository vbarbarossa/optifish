

l1 = 20
l2 = 10
inter_area <- c(l1,l2)
passability_dam <- c(1,0.5)

m = matrix(1,nrow = 2, ncol = 2)
m = m*inter_area
m = t(m)*inter_area
m = m/sum(inter_area)**2

# passability matrix
pm = matrix(1,nrow = 2, ncol = 2)

# for(i in 1:ncol(pm)) for(j in 1:nrow(pm)) pm[j,i] <- pm[j,i]*passability_dam[j]
# m = t(m)*inter_area
# m = m/sum(inter_area)**2

HYBAS <- data.frame(INTER_ID = c(1,2), INTER_NEXT = c(0,1))
df <- data.frame(from = HYBAS$INTER_ID,to = HYBAS$INTER_NEXT)
outlet <- as.character(df$from[which(df$to == 0)])

df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g <- graph_from_data_frame(d=df)
I <- get.adjacency(g)

diag(I) <- 1

for(i in 1:nrow(I)) I[i,1] <- I[i,1] * passability_dam[i]
I = I + t(I)
diag(I) <- 1

l1 = 15
l2 = 10
l3 = 5

C = matrix(c(1,0.5,0.3,0.5,1,0.15,0.3,0.15,1),nrow=3)
I = matrix(1, nrow = 3, ncol = 3)
p = c(l1,l2,l3)
c = c(1,0.5,0.3)

L = (I*p * t(I*p)) / sum(p)**2 
C = (I*c * t(I*c))
diag(C) <- 1

sum(C*L)

C1 <- C
C1[1,3] <- C1[3,1] <- 1

sum(C1*L)

C1[2,3] <- C1[3,2] <- 0.5



library(igraph)
# try with similar structure data
HYBAS <- data.frame(INTER_ID = c(1,2,3,4), INTER_NEXT = c(0,1,1,3))
df <- data.frame(from = HYBAS$INTER_ID,to = HYBAS$INTER_NEXT)
outlet <- as.character(df$from[which(df$to == 0)])

df[which(df$to == 0),]$to <- outlet %>% as.numeric()
g <- graph_from_data_frame(d=df)
plot(g)
E(g)$pass <- c(1,0.5,0.3,0.3)

vertices_id <- names(igraph::V(g))
C <- 10^dodgr::dodgr_dists(
  mutate(
    rbind( 
      igraph::as_data_frame(g, what = "edges"), 
      igraph::as_data_frame(as.undirected(g,mode = 'each'), what = "edges")
    ),dist = log10(.data$pass)
  ),
  from = vertices_id, to = vertices_id)

l = c(15,10,3,2)
I = matrix(1, nrow = length(l), ncol = length(l))
L = (I*l * t(I*l)) / sum(l)**2 

sum(C*L)
