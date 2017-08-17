# Name: vectors.R
# Auth: uhkniazi
# Date: 09/08/2017
# Desc: vectors, matrices and others - blog related script

p.old = par()

df = data.frame(x = c(-3,-3,-1,-1,0,0,1,2,2,3), 
                y = c(-8,-5,-3,0,-1,0,5,1,6,5))
row.names(df) = paste0('S', 1:10)

par(mfrow=(c(1,2)))

plot(df$x, df$y, xlim=c(-10, 10), ylim=c(-10, 10), main='Centered Data in Variable Space',
     xlab='variable 1 (x)', ylab='variable 2(y)', pch=19)

df2 = data.frame(t(df))

plot(df2$S1, df2$S2, main='First 2 data points in Subject Space',
     xlab='Subject 1 (S1)', ylab='Subject 2(S2)', pch=19)
text(df2$S1, df2$S2, labels = c('x', 'y'), adj = 2, pos = c(1,3))

rad2deg = function(rad) {(rad * 180) / (pi)}
deg2rad = function(deg) {(deg * pi) / (180)}

cor(df$x, df$y); angleRad = acos(cor(df$x, df$y)); 
rad2deg(acos(cor(df$x, df$y)))

## coordinates for vector X
magX = norm(df$x, '2')
ptX = c(x=magX * cos(0), y=magX * sin(0))

magY = norm(df$y, '2')
ptY = c(x=magY * cos(angleRad), y=magY * sin(angleRad))

apply(df, 2, sd)

lincomb = function(x, y){
  return(3*x - 2*y)
}

lincomb(c(2,1), c(-1,3))

## define 2 vectors
x = c(2,1); y = c(-1, 3)

## dot product 
t(x) %*% y
## magnitude x
norm(x, '2')
## or
sqrt(t(x) %*% x)
## magnitude y
norm(y, '2')
## calculate cosine of angle
1/( 2.236068 * 3.162278)
## convert this to the angle
acos(0.1414213)
## convert it to degrees
rad2deg(1.428899)

### input output operations, functions, vectors and matrices
## some operations or functions
op1 = function(x, y, z){
  return(3*x + 4*y + 5*z)
}

op2 = function(x, y, z){
  return(3*x + 0*y + 0*z)
}

op3 = function(x, y, z){
  return(4*x + (-3)*y + 2*z)
}


## some vector inputs
inp1 = c(5, 10, 12)
inp2 = c(4, 6, -2)
inp3 = c(8, -5, 4)

## perform operation on input
op1(inp1[1], inp1[2], inp1[3])
op2(inp2[1], inp2[2], inp2[3])
op3(inp3[1], inp3[2], inp3[3])

## create input matrix
mInputs = cbind(inp1, inp2, inp3)

## operations matrix
mOperations = rbind(c(3, 4, 5),
                    c(3, 0, 0),
                    c(4, -3, 2))

## matrix multiplication
## operations %*% inputs
mOutput = mOperations %*% mInputs
## resultant matrix = rows=size of operation (2) and cols=number of inputs (2)
mOutput

## what is the determinant of the matrix
det(mOperations)
## get inverse of the matrix
mOperations.inv = solve(mOperations)
## retrieve original input
mOperations.inv %*% mOutput









