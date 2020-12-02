import File
FileName = "Data.txt"
WirteData = ["NJ"]

def writeMatrix(matrix,length,WirteData):
    for i in range(length):
        ans = ""
        for j in range(length):
            ans = ans+str(matrix[(length)*(i)+j])+","
        ans=ans[0:-1]
        WirteData.append(ans)

def DisplayMatrix(matrix,col,row):
    for i in range(row):
        for j in range(col):
            print(matrix[(col)*(i)+j],end=" ")
        print()

def D2Tri(matrix,length):
    trans = []
    for i in range(1,length):
        for j in range(i,1+length):
            trans.append(-1) 
        for j in range(1,i+1):
            trans.append(matrix[length*i+j-1])
    trans.append(-1)
    return(trans)

def mincell(matrix,length):
    ans = (matrix[length*1+1-1])
    loc = length*1+1-1    
    for i in range(1,length):
        for j in range(1,i+1):
            if (ans > (matrix[length*i+j-1])):
                ans = (matrix[length*i+j-1])
                loc = length*i+j-1
    return([ans,loc])

def M2Q(matrix,length,trans):
    Q=[-1]*length*length
    for Row in range(length):
        for Col in range(length):
            if trans[Row*length+Col] != -1:
                ans = (length-2)*matrix[Row*length+Col]
                for m in range(length):
                    ans = ans - matrix[Col*length+m] - matrix[Row*length+m]
                Q[Row*length+Col] = ans 
    return(Q)
def Antitrans(trans,length):
    anti =  trans 
    for Row in range(length):
        for Col in range(length):
            if (Row == Col):
                anti[Row*length+Col] = 0
            elif (Col>Row):
                anti[Row*length+Col] = anti[Col*length+Row] 
    return anti

def onestep(matrix,trans,Q,length,name,namenum,distance,distancename):
    WirteData.append("n = "+str(length))
    WirteData.append("Distance Matrix")
    WirteData.append(str(name))
    writeMatrix(matrix,length,WirteData)
    WirteData.append("\n")
    WirteData.append("Q Matrix")
    WirteData.append(str(name))
    writeMatrix(Q,length,WirteData)
    WirteData.append("\n")
    tem = matrix
    matrix = Q
    Q = tem
    [mindata,minloc] = mincell(matrix,length)
    Row = ((int)(minloc/length))
    Col = ((int)(minloc%length))
    if length != 2:
        ans = 0
        ans2 = 0.5*Q[Row*length+Col]
        for m in range(length):
            Rowtem2 = max(m,Col)
            Coltem2 = min(m,Col)
            Rowtem3 = max(m,Row)
            Coltem3 = min(m,Row)
            ans = ans + (trans[Rowtem2*length+Coltem2] - trans[Rowtem3*length+Coltem3])/(2*(length-2))
        delta = (ans2 - ans)
        if (delta>Q[Row*length+Col]):
            delta = ans2+ans
    else:
        ans2 = Q[Row*length+Col]
        
    NewMatrix = []
    for i in range(len(matrix)):
        if ((int)(i/length) != Row and (int)(i%length) != Col and (int)(i/length) != Col and (int)(i%length) != Row):
            # print("No")
            NewMatrix.append(Q[i])

    New = []
    for i in range(length):
        if (i != Col and i != Row):
            Rowtem = max(i,Row)
            Coltem = min(i,Row)
            d1 = (Q[Rowtem*length+Coltem])
            # print(d1)
            Rowtem = max(i,Col)
            Coltem = min(i,Col)
            d2 = (Q[Rowtem*length+Coltem])
            # print(d2)
            Rowtem = max(Row,Col)
            Coltem = min(Col,Col)
            d3 = (Q[Rowtem*length+Coltem])
            NewEle = (d1+d2-d3)/2
            NewMatrix.append(NewEle)
    length = length - 1
    for i in range(length):
        NewMatrix.insert((i+1)*length-1,-1)
    temmax = max(Row,Col)
    temmin = min(Row,Col)
    if length != 1:
        namenum.append(namenum[Row]+namenum[Col])
        name.append([name[Row],name[Col]])
        distance.append([ans2-ans,ans2+ans])
        distancename.append([name[Row],name[Col]])
        name.pop(temmax)
        name.pop(temmin)
        namenum.pop(temmax)
        namenum.pop(temmin)
        matrix = NewMatrix
        trans = Antitrans(matrix,length)
        Q = (M2Q(trans,length,matrix))
        Q = D2Tri(Q,length)
    if length != 1:
       [distance,distancename]=onestep(matrix,trans,Q,length,name,namenum,distance,distancename)
    else:
        namenum.append(namenum[Row]+namenum[Col])
        name.append([name[Row],name[Col]])
        distance[-1].append(ans2)
        distancename[-1].append(name[Col])
        name.pop(temmax)
        name.pop(temmin)
        namenum.pop(temmax)
        namenum.pop(temmin)
    return(distance,distancename)

data = File.OpenFile(FileName)
name = []
length = len(data)-1
namenum = [1]*length
matrix = [0]*(length*length)
for i in range(length+1):
    if i == 0:
        name = data[i][0:-1].split('\t')
    else:
        line = data[i][0:-1].split('\t')
        for j in range(length):
            matrix[(i-1)*length+j] = float(line[j])
# print(matrix)
distance = []
distancename = []
trans = D2Tri(matrix,length)
Q = (M2Q(matrix,length,trans))
tem = matrix 
matrix =  D2Tri(matrix,length)

[distance,distancename] = onestep(matrix,tem,Q,length,name,namenum,distance,distancename)
final = distancename[-1]
tem = str(final)
for i in range(len(distancename)-1,-1,-1):
    loc = (tem.find(str(distancename[i])))
    if(len(distance[i])==3):
        rep = "("+str(distancename[i][0])+":"+str(distance[i][0])+","+str(distancename[i][1])+":"+str(distance[i][1])+","+str(distancename[i][2])+":"+str(distance[i][2])+")"
    else:
        rep = "("+str(distancename[i][0])+":"+str(distance[i][0])+","+str(distancename[i][1])+":"+str(distance[i][1])+")"
    tem = tem[0:loc]+rep+tem[loc+len(str(distancename[i]))::]

WirteData.append("\n")
WirteData.append("Final Newick Format")
WirteData.append(str(tem))
File.WriteFile("NJ.txt",WirteData)