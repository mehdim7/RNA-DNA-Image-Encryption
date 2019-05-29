'''
Image Encryption
Combination RNA and DNA and Logistic map
CopyRight Mehdi Yadollahi  January 5, 2019'''
import random
import cv2
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from numpy import unique
from scipy.stats import entropy as scipy_entropy
import time

class Gen:
    def __init__(self,bases,code):
        self.bases=bases
        self.code=code

    def getCode(self):
        return self.code

    def getBases(self):
        return self.bases

class Codon:
    def __init__(self,gen1,gen2,gen3):
        self.codon=[gen1,gen2,gen3]

    def getString(self):
        s=''
        for i in range(3):
            s=s+self.codon[i].getBases()
        return s

    def getBinary(self):
        s = ''
        for i in range(3):
            s = s + self.codon[i].getCode()
        return s

    def getComplement(self,basesTemp):
        newCodon=[]
        for i in range(3):
            if(self.codon[i].bases=='A'):
                newCodon.append(basesTemp[3])
            elif (self.codon[i].bases == 'U'):
                newCodon.append(basesTemp[0])
            elif (self.codon[i].bases == 'C'):
                newCodon.append(basesTemp[2])
            elif (self.codon[i].bases == 'G'):
                newCodon.append(basesTemp[1])
            else:
                print('Error complementary')
        return Codon(newCodon[0],newCodon[1],newCodon[2])

    def Xor(self,temp):
      #  print('>',temp.getBinary())
       # print('^',self.getBinary())
        bin3=int(self.getBinary(),2) ^ int(temp.getBinary(),2)
        bin3='{0:06b}'.format(bin3)
        if len(bin3)!=6:
            print('XOR Error')
     #   print('=',bin3)

        return bin3


class CodonRow:
    def __init__(self,indexC,randC,codonC):
        self.index=indexC
        self.rand=randC
        self.codons=codonC

    def getRow(self):
        return [self.index,self.rand,self.codons]

    def getCodon(self):
        return self.codons

    def getIndex(self):
        return  self.index

    def getRand(self):
        return self.rand

class DNA:
    def __init__(self,imgBit,key):
        self.DNATableBits={'A':['00','00','11','11','10','01','10','01'],
               'T': ['11', '11', '00', '00', '01', '10', '01', '10'],
               'C': ['10','01','10','01', '00','00','11','11'],
               'G': ['01', '10', '01','10','11', '11', '00', '00']}

        self.DNATable = [ {'00':'A', '11':'T', '10':'C', '01':'G'},
                           {'00': 'A', '11': 'T', '01': 'C', '10': 'G'},
                           {'11': 'A', '00': 'T', '10': 'C', '01': 'G'},  #2
                           {'11': 'A', '00': 'T', '01': 'C', '10': 'G'},
                           {'10': 'A', '01': 'T', '00': 'C', '11': 'G'},  #4
                           {'01': 'A', '10': 'T', '00': 'C', '11': 'G'},
                           {'10': 'A', '01': 'T', '11': 'C', '00': 'G'},  #6
                           {'01': 'A', '10': 'T', '11': 'C', '00': 'G'}]

        self.XORDNATable={'A':{'A':'A' , 'T':'T' , 'C':'C', 'G':'G'},
                          'T':{'A':'T','T':'A','C':'G','G':'C'},
                          'C':{'A':'C','T':'G','C':'A','G':'T'},
                          'G':{'A':'G', 'T':'C' , 'C':'T' , 'G':'A'}}

        self.imgBit=imgBit
        self.key=key
        self.R=3.999
        self.x0=self.createKey(key)
        print('Key: ',self.key, 'X0: ',self.x0)


    def createKey(self,key):
        a = ord(key[0])
        for i in range(1, 12):
            a = a ^ ord(key[i])
        a = a / (2 ** 12)
        b = ''
        for i in range(12, len(key)):
            b += '{0:08b}'.format(ord(key[i]))
        b = int(b, 2) / (2 ** 32)
        return (a + b) / 2

    def DNACrypt(self):
        finalBit = ''
        keyCounter = 0
        lastPixcel='00000000'
        logistic=self.logisticMap()
        for i in range(0, len(self.imgBit), 8):
            x1 = round(next(logistic)*255) #random.randint(0, 255)   #key
            x2 = round(next(logistic)*7)
            x3 = round(next(logistic)*7)
            x4 = round(next(logistic)*7)
            x5 = round(next(logistic)*7)

            key8Bit='{0:08b}'.format(x1)
            data8Bit=self.imgBit[i:i + 8]
            key4DNA=self.get4DNA(key8Bit,x2)   # key => x2
            data4DNA=self.get4DNA(data8Bit,x3)   # data => x3
            result4DNA=self.XorDNA(key4DNA,data4DNA)
            lastPixcel4DNA=self.get4DNA(lastPixcel,x4)    # lastPixcel => x4
            result4DNA=self.XorDNA(result4DNA,lastPixcel4DNA)

            lastPixcel=self.get4DNABinary(result4DNA,x5)    # result => x5
            finalBit += lastPixcel

        return finalBit

    def logisticMap(self):
        x=self.x0
        while True:
            x = self.R * x * (1.0 - x)
            yield x

    def get4DNABinary(self,result4DNA,ruleNumber):
        cc = ''
        for i in range(0, len(result4DNA)):
            cc += self.DNATableBits[result4DNA[i]][ruleNumber]
        return cc

    def get4DNA(self,bits,ruleNumber):
        cc=''
        for i in range(0, len(bits), 2):
            temp = bits[i:i + 2]
            cc +=self.DNATable[ruleNumber][temp]
        return cc

    def XorDNA(self,result4DNA,lastPixcel4DNA):
        cc=''
        for i in range(0, len(result4DNA)):
            cc+=self.XORDNATable[result4DNA[i]][lastPixcel4DNA[i]]
        return cc

class RNA:
    baseRNACodonTable=[]
    def __init__(self,imgBit,keyBit):
        self.imgBit=imgBit
        self.key=keyBit
       # self.key='00110110'
        self.gen1=  Gen('A','00')
        self.gen2 = Gen('C', '01')
        self.gen3 = Gen('G', '10')
        self.gen4 = Gen('U', '11')
        self.Bases=[self.gen1,self.gen2,self.gen3,self.gen4]
        self.baseRNACodonTable= self.CreateRNACodonTable()
        self.RNACodonTable1= self.CreateRNACodonTable(setRand=True)
        self.RNACodonTableC1= self.CreateRNACodonTableC(self.RNACodonTable1)
        self.RNACodonTable2= self.CreateRNACodonTable(setRand=True)
        self.RNACodonTableC2= self.CreateRNACodonTableC(self.RNACodonTable2)

       # self.ShowCodonTable(self.baseRNACodonTable)

    def CreateRNACodonTable(self,setRand=False):
        index=0
        rand = 1
        if(setRand==True):
            rand=random.uniform(0, 1)

        codonTable=[]
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    newRow=CodonRow(index,rand,Codon(self.Bases[i], self.Bases[j], self.Bases[k]))
                    codonTable.append(newRow)
                    index+=1
                    rand=index
                    if (setRand == True):
                        rand=random.uniform(0, 1)

        codonTable.sort(key=lambda x: x.rand, reverse=False)
        return codonTable

    def CreateRNACodonTableC(self, codonTable):
        codonTableC=[]
        for i in range(64):
            tt= codonTable[i].getCodon().getComplement(self.Bases)
            newRow=CodonRow(codonTable[i].index,codonTable[i].rand,tt)
            codonTableC.append(newRow)
        return codonTableC

    def ShowCodonTable(self,codonTable):
        for i in range(64):
            st=str.format("{0},{1},{2},{3},{4}",i,codonTable[i].getIndex(),codonTable[i].getRand()
                  ,codonTable[i].getCodon().getString()
                  ,codonTable[i].getCodon().getBinary())
            print(st)

    def getRNACodeCodon(self,tripleCodon,key2Bit):
        if(key2Bit=='00'):
            table=self.RNACodonTable1
        elif(key2Bit=='01'):
            table=self.RNACodonTableC1
        elif (key2Bit == '10'):
            table = self.RNACodonTable2
        elif (key2Bit == '11'):
            table = self.RNACodonTableC2
        else:
            print('Error getRNACodeCodon:',key2Bit,tripleCodon)
        indexTemp=self.findTripleCodon(tripleCodon)
        return indexTemp,table[indexTemp].getCodon()

    def findTripleCodon(self,tripleCodon):
        for i in range(64):
            if (tripleCodon==self.baseRNACodonTable[i].getCodon().getString()):
                return i
        print(tripleCodon)
        print('Error finding triple Codon')
        return -1

    def convertBinaryToBases(self,bits):
        cc = ''
        for i in range(0, len(bits), 2):
            temp = bits[i:i + 2]
            cc += list(filter(lambda x: x.getCode() == temp, self.Bases))[0].getBases()

        fix = len(cc) % 3
        fix= (3-fix)% 3
        p = 'A' * fix
        cc += p
        return cc

    def getXorCodon(self,tempCodon,sourceCodon):
        bin1=tempCodon.Xor(sourceCodon)
        tripleCodon=self.convertBinaryToBases(bin1)
        tempindex= self.findTripleCodon(tripleCodon)
        return self.baseRNACodonTable[tempindex].getCodon()

    def RNACrypt(self):
        sourceString=self.convertBinaryToBases(self.imgBit)
        finalString=''
        finalBit=''

        keyCounter=0
        for i in range(0,len(sourceString),3):
            startKey=keyCounter%(len(self.key)-1)
            keyCounter+=1
            tt=self.key[startKey:startKey+2]
            index,tempCodon= self.getRNACodeCodon(sourceString[i:i + 3], tt)
            finalString+=tempCodon.getString()
            finalBit+=tempCodon.getBinary()

        return finalBit

class ImageCrypt:
    def __init__(self,img,pathRNAFinal,key):
        self.imgSource=img
        self.pathRNAFinal=pathRNAFinal

        self.imgRNAFinal=img.copy()

        finalBit=self.convertImageToBinary(self.imgSource)

        self.dna=DNA(finalBit,key)
        finalBit = self.dna.DNACrypt()

        self.rna = RNA(finalBit, self.convertKeyToBinary(key))
        finalBit = self.rna.RNACrypt()

        self.saveRNAFinalImage(finalBit)

    def convertKeyToBinary(self,key):
        a = ''
        for i in range(0, len(key)):
            a += '{0:08b}'.format(ord(key[i]))
        return a

    def convertImageToBinary(self,img):
        bb = ''
        rows, cols = img.shape
        for i in range(rows):
            for j in range(cols):
                bb += '{0:08b}'.format(img[i, j])
        return bb

    def saveRNAFinalImage(self,finalBit):
        rows, cols = self.imgRNAFinal.shape
        for i in range(rows):
            for j in range(cols):
                start=i * cols * 8 + j * 8
                temp8Bit= finalBit[start:start+8]
                self.imgRNAFinal[i,j] = int(temp8Bit,2)
        cv2.imwrite(self.pathRNAFinal, self.imgRNAFinal)

class DataSet:
    def __init__(self,imagesDict):
        self.images=imagesDict
        self.path=path
        self.filePass='.bmp'

    def getDataSetPath(self,dataSetName,filepass='.jpg'):
        return self.path+self.images[dataSetName]+filepass

    def getDataSetFinalPath(self,dataSetName):
        return self.path + self.images[dataSetName] +'RNA'+ self.filePass

    def getDataSetFinalPathKeyChange(self,dataSetName):
        return self.path + self.images[dataSetName] +'RNAKeyChange'+ self.filePass

    def getDataSetFinalPathKeyChangeResult(self,dataSetName,count):
        return self.path + self.images[dataSetName] +'KeyChange_'+str(count)+ self.filePass

    def getDataSetFinalPathUACI1(self,dataSetName):
        return self.path + self.images[dataSetName] +'RNAUACI1'+ self.filePass

    def getDataSetFinalPathUACI2(self,dataSetName):
        return self.path + self.images[dataSetName] +'RNAUACI2'+ self.filePass

    def getDataSetFinalPathUACI3(self,dataSetName):
        return self.path + self.images[dataSetName] +'RNAUACI3'+ self.filePass

    def getDataSetRNAPlt(self,dataSetName,filepass='.png'):
        return self.path + self.images[dataSetName] +'RNAPlt'+ filepass

    def getDataSetHorizontalCorPlt(self,dataSetName,filepass='.png'):
        return self.path + self.images[dataSetName] +'HorizontalCorPlt'+ filepass

    def getDataSetVerticalCorPlt(self,dataSetName,filepass='.png'):
        return self.path + self.images[dataSetName] +'VerticalCorPlt'+ filepass

    def getDataSetDiagonalCorPlt(self, dataSetName, filepass='.png'):
        return self.path + self.images[dataSetName] + 'DiagonalCorPlt' + filepass


class AnalyzeImage:

    def __init__(self,dataSetName,key,dataSets,infoFile):
        self.dataSets=dataSets
        self.dataSetName=dataSetName
        self.infoFile=infoFile
        self.key=key
        self.correlationcCount=8000

        self.imgSource = cv2.imread(self.dataSets.getDataSetPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        self.imgRNAFinal = cv2.imread(self.dataSets.getDataSetFinalPath(dataSetName), cv2.IMREAD_GRAYSCALE)

        entropyFinal = self.Entropy(self.imgRNAFinal)
        print('RNA Entropy: ', entropyFinal)
        self.infoFile.write('\nRNA Entropy: ' + str(entropyFinal))

        self.key2 = self.getNextKey(self.key)
        run = ImageCrypt(self.imgSource, self.dataSets.getDataSetFinalPathKeyChange(dataSetName), self.key2)

        count, countPercent =self.keySensitivityImage(dataSetName)
        infoFile.write('\nKey Sensivity: ' + str(count)+' , '+ str(round(countPercent,2)))
        self.NPCRUACI(dataSetName, infoFile)
        self.Correlation()


    def Entropy(self,img):
        _, counts = unique(img, return_counts=True)
        return scipy_entropy(counts, base=2)

    def Correlation(self):
        matrixX,matrixY=self.getCorrelation(self.imgSource, 1, 0)
        matrixX2,matrixY2=self.getCorrelation(self.imgRNAFinal, 1, 0)
        self.infoFile.write('\nCorrelation Coefficient:')
        corHorizontal = np.corrcoef(matrixX2, matrixY2)
        self.infoFile.write('\nHorizontal: ' + str(round(corHorizontal[1, 0],4)))

        matplotlib.style.use('ggplot')
        plt.subplot(122), plt.scatter(matrixX2, matrixY2,c='blue', alpha=0.5), plt.title('Horizontal')
        plt.subplot(121), plt.scatter(matrixX, matrixY,c='blue', alpha=0.5), plt.title('ORIGINAL')
        plt.savefig(self.dataSets.getDataSetHorizontalCorPlt(self.dataSetName))
        #plt.show()

        matrixX, matrixY = self.getCorrelation(self.imgSource, 0, 1)
        matrixX2, matrixY2 = self.getCorrelation(self.imgRNAFinal, 0, 1)

        corVertical = np.corrcoef(matrixX2, matrixY2)
        self.infoFile.write('\nVertical: ' + str(round(corVertical[1, 0],4)))

        matplotlib.style.use('ggplot')
        plt.subplot(122), plt.scatter(matrixX2, matrixY2,c='blue', alpha=0.5), plt.title('Vertical')
        plt.subplot(121), plt.scatter(matrixX, matrixY,c='blue', alpha=0.5), plt.title('ORIGINAL')
        plt.savefig(self.dataSets.getDataSetVerticalCorPlt(self.dataSetName))

        matrixX, matrixY = self.getCorrelation(self.imgSource, 1, 1)
        matrixX2, matrixY2 = self.getCorrelation(self.imgRNAFinal, 1, 1)

        corDiagonal = np.corrcoef(matrixX2, matrixY2)
        self.infoFile.write('\nDiagonal: ' + str(round(corDiagonal[1,0],4)))

        matplotlib.style.use('ggplot')
        plt.subplot(122), plt.scatter(matrixX2, matrixY2,c='blue', alpha=0.5), plt.title('Diagonal')
        plt.subplot(121), plt.scatter(matrixX, matrixY,c='blue', alpha=0.5), plt.title('ORIGINAL')
        plt.savefig(self.dataSets.getDataSetDiagonalCorPlt(self.dataSetName),dpi=300)

    def getCorrelation(self,img,x,y):   # x=1 y=0  horizontal    x=0 y=1 vertical    x=1 y=1 diagonal
        matrixX=[]
        matrixY=[]

        rows, cols = img.shape

        maxRow=rows-1
        maxCol=cols-1
        if x==0 and y==1:
            maxCol=cols
        if x==1 and y==0:
            maxRow=rows

        random.seed(1)
        for i in range(self.correlationcCount):
            rx=random.randint(0,maxRow-2)
            ry=random.randint(0,maxCol-2)
            matrixX.append(img[rx,ry])
            matrixY.append(img[rx+x,ry+y])

        return matrixX,matrixY

    def getNextKey(self, key):
        if (len(key) < 16):
            print('Error key')
        if (ord(key[15]) == 255):
            key = key[0:15] + chr(0)
        else:
            key = key[0:15] + chr(ord(key[15]) + 1)
        return key

    def getNextPixcel(self, pixcel):
        if (pixcel == 255):
            pixcel = 0
        else:
            pixcel += 1
        return pixcel

    def keySensitivityImage(self, dataSetName):
        stt=self.dataSets.getDataSetFinalPath(dataSetName)
        img1 = cv2.imread(self.dataSets.getDataSetFinalPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        img2 = cv2.imread(self.dataSets.getDataSetFinalPathKeyChange(dataSetName), cv2.IMREAD_GRAYSCALE)
        sensitivityImage = cv2.imread(self.dataSets.getDataSetFinalPathKeyChange(dataSetName), cv2.IMREAD_GRAYSCALE)
        rows, cols = img1.shape
        count = 0
        for i in range(rows):
            for j in range(cols):
                if (img1[i, j] == img2[i, j]):
                    sensitivityImage[i, j] = 255
                    count += 1
                else:
                    sensitivityImage[i, j] = 0

        cv2.imwrite(self.dataSets.getDataSetFinalPathKeyChangeResult(dataSetName, count), sensitivityImage)
        return count,((1-(count/(rows*cols)))*100)

    def NPCRUACI(self, dataSetName, infoFile):
        img1 = cv2.imread(self.dataSets.getDataSetPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        img1[0, 0] = self.getNextPixcel(img1[0, 0])
        rows, cols = img1.shape
        run1 = ImageCrypt(img1, self.dataSets.getDataSetFinalPathUACI1(dataSetName), self.key)

        NPCR1, UACI1 = self.getEachNPCRUACI(self.dataSets.getDataSetFinalPath(dataSetName),
                                        self.dataSets.getDataSetFinalPathUACI1(dataSetName))
        infoFile.write('\nFirst Pixcel: ')
        infoFile.write('\nNPCR: ' + str(round(NPCR1,6)))
        infoFile.write('\nUACI: ' + str(round(UACI1,6)))

        img2 = cv2.imread(self.dataSets.getDataSetPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        img2[int(rows / 2), int(cols / 2)] = self.getNextPixcel(img2[int(rows / 2), int(cols / 2)])
        run2 = ImageCrypt(img2, self.dataSets.getDataSetFinalPathUACI2(dataSetName), self.key)

        NPCR2, UACI2 = self.getEachNPCRUACI(self.dataSets.getDataSetFinalPath(dataSetName),
                                        self.dataSets.getDataSetFinalPathUACI2(dataSetName))
        infoFile.write('\nMidle Pixcel: ')
        infoFile.write('\nNPCR: ' + str(round(NPCR2,6)))
        infoFile.write('\nUACI: ' + str(round(UACI2,6)))

        img3 = cv2.imread(self.dataSets.getDataSetPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        img3[rows - 1, cols - 1] = self.getNextPixcel(img2[rows - 1, cols - 1])
        run3 = ImageCrypt(img3, self.dataSets.getDataSetFinalPathUACI3(dataSetName), self.key)

        NPCR3, UACI3 = self.getEachNPCRUACI(self.dataSets.getDataSetFinalPath(dataSetName),
                                        self.dataSets.getDataSetFinalPathUACI3(dataSetName))
        infoFile.write('\nLast Pixcel: ')
        infoFile.write('\nNPCR: ' + str(round(NPCR3,6)))
        infoFile.write('\nUACI: ' + str(round(UACI3,6)))

        NPCR4 = (NPCR1 + NPCR2 + NPCR3) / 3
        UACI4 = (UACI1 + UACI2 + UACI3) / 3

        infoFile.write('\nAverage: ')
        infoFile.write('\nNPCR: ' + str(round(NPCR4,6)))
        infoFile.write('\nUACI: ' + str(round(UACI4,6)))

    def getEachNPCRUACI(self, imgPath1, imgPath2):
        img1 = cv2.imread(imgPath1, cv2.IMREAD_GRAYSCALE)
        img2 = cv2.imread(imgPath2, cv2.IMREAD_GRAYSCALE)

        rows, cols = img1.shape
        sumNPCR = 0
        sumUACI = 0
        for i in range(rows):
            for j in range(cols):
                sumUACI += abs(int(img1[i, j]) - int(img2[i, j]))
                if (img1[i, j] != img2[i, j]):
                    sumNPCR += 1
        NPCR = sumNPCR / (rows * cols)
        UACI = sumUACI / (255 * rows * cols)
        return NPCR, UACI


class ProposedAlgorithm:
    def __init__(self,dataSets,key):
        self.dataSets=DataSet(dataSets)
        self.key=key
        for k in dataSets:
            print('>>>> ', k, '  <<<<<')
            start=time.time()
            infoFile = open(self.dataSets.getDataSetPath(k,'.txt'),'w')
            infoFile.write('Image: '+k)
            #infoFile.write('\r\n'+'Key: ' + self.key)
            img1 = cv2.imread(self.dataSets.getDataSetPath(k), cv2.IMREAD_GRAYSCALE)
            run=ImageCrypt(img1,self.dataSets.getDataSetFinalPath(k),self.key)
            analyze=AnalyzeImage(k,self.key,self.dataSets,infoFile)
            #self.showImage(k)
            end=time.time()
            infoFile.write('\nTime: ' + str(round(end-start,3)))
            infoFile.close()

    def showImage(self,dataSetName):
        imgSource =cv2.imread(self.dataSets.getDataSetPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        imgRNA=cv2.imread(self.dataSets.getDataSetFinalPath(dataSetName), cv2.IMREAD_GRAYSCALE)
        imgKeySensitive=cv2.imread(self.dataSets.getDataSetFinalPathKeyChange(dataSetName), cv2.IMREAD_GRAYSCALE)

        plt.subplot(231), plt.imshow(imgSource), plt.title('ORIGINAL')
        plt.subplot(232), plt.imshow(imgKeySensitive), plt.title('Key Sensitive')
        plt.subplot(233), plt.imshow(imgRNA), plt.title('RNA')
        hist1 = cv2.calcHist([imgSource], [0], None, [256], [0, 256])
        hist2 = cv2.calcHist([imgKeySensitive], [0], None, [256], [0, 256])
        hist3 = cv2.calcHist([imgRNA], [0], None, [256], [0, 256])

        plt.subplot(234), plt.hist(imgSource.ravel(), 256, [0, 256])
        plt.subplot(235), plt.hist(imgKeySensitive.ravel(), 256, [0, 256])
        plt.subplot(236), plt.hist(imgRNA.ravel(), 256, [0, 256])

        plt.savefig(self.dataSets.getDataSetRNAPlt(dataSetName))


# key for image encryption
key='Ì¼ÆN²ÜYTö,¡Ä'          # 128 bit ( 16 char)
# Name of DataSets
dataSets={'10':'10','11':'11','12':'12','13':'13','14':'14','15':'15'
    ,'16': '16','17':'17','18':'18','19':'19'}


# DataSet Directory Path
path='D:/RNA/'
p=ProposedAlgorithm(dataSets,key)
print('Done!!!')
