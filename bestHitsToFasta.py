#!/usr/bin/python3.6.4
# -*- coding: utf-8 -*-
### input : output de diamond outfmt -6
### sort according to decrease bitscore 
### ouput : Best hits no overlap, no duplicats
### collapsing of near fragments less than 50 nucl

##--- Header
import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import sys      # gestion arguments et opt

inFileNoSort= sys.argv[1]
inFile= sys.argv[1]+"_sorted"
outFile=sys.argv[1]+".bed1"


##--- Main
os.system("sort -rgk'12,12' "+inFileNoSort+' > '+inFile ) ## decreasing sort

# Open *lecture*:
blast = open(inFile, "r")   # read input
out = open(outFile, "w")    # write output
bestHits={}     # dico of hits


for ligne in blast: # --- Parsing of input ( sorting blast )
    ligne=ligne.rstrip()    # delete sauts de lignes
    l=ligne.split("\t")     # split line on tab into list
    if (int(l[6]) < int(l[7])):   ## test sens du hit            
        deb = l[6]
        fin = l[7]
        stranded='+'
    else:               # si hit sur le brin revers
        deb = l[7]
        fin = l[6]
        stranded='-'    
   # stock les eve dans un dico (key=nom de scaffold/value = liste des eve du scaffold)
    if l[0] in bestHits: ## if scaffold in dico
        ll=0     # flag : pas modif avec un eve deja existant dans le dico
        for i in range (0,len(bestHits[l[0]])): # pour chaque eve du scaffold present dans dico
            if (int(deb) >= int(bestHits[l[0]][i]['start'])) and (int(fin) <= int(bestHits[l[0]][i]['stop'])) : # si eve_Ligne inclu eve_dico
                ll=ll+1
                break     #  modif flag = inclu (eve_Ligne ds eve_dico)
            elif (int(deb) <= int(bestHits[l[0]][i]['start'])) and (int(fin) >= int(bestHits[l[0]][i]['stop'])) : # si il est inclu
                bestHits[l[0]][i]['start']= deb    #  modif start et stop = eve_dico inclu dans eve_Ligne
                bestHits[l[0]][i]['stop']= fin
                ll+=1
                break
            elif (int(deb) < int(bestHits[l[0]][i]['start'])) and (int(fin) >= int(bestHits[l[0]][i]['start'])) and (int(fin) <= int(bestHits[l[0]][i]['stop'])): # si overlap en debut du hits
                bestHits[l[0]][i]['start']= deb     #  modif start = olap fin eve_Ligne
                ll+=1
                break
            elif (int(fin) > int(bestHits[l[0]][i]['stop'])) and (int(deb) >= int(bestHits[l[0]][i]['start'])) and (int(deb) <= int(bestHits[l[0]][i]['stop'])): # si overlap a la fin du hits
                bestHits[l[0]][i]['stop']= fin     #  modif stop = olap debut eve_Ligne
                ll=ll+1
                break
            else:   ## ni inclu, ni overlap
                distStart= int(bestHits[l[0]][i]['start']) - int(fin)   # distance si match avant hit
                distStop = int(deb) - int(bestHits[l[0]][i]['stop'])     # distance si match apres hit
                if ((distStart <= 100)  and (distStart > 0 )  and  (deb < bestHits[l[0]][i]['start'])) :
                    bestHits[l[0]][i]['start']= deb # modif debut de l'eve
                    ll=ll+1
                    break
                  
                elif ((distStop <= 100)  and (distStop > 0) and (fin > bestHits[l[0]][i]['stop'])) :
                    bestHits[l[0]][i]['stop']= fin # modif fin de l'eve
                    ll=ll+1
                    break 
        if (ll ==0) : ## si aucune modif on ajout la region    
            bestHits[l[0]].append({'contig':l[0],'start':deb,'stop':fin,'subj':l[1],'subj-start':l[8],'subj-stop':l[9],'id':l[2],'Evalue':l[10],'bitScore':l[11],'len':l[3],'sens':stranded })
    else : ## nouveau scaffold
        bestHits[l[0]] = [{'contig':l[0],'start':deb,'stop':fin,'subj':l[1],'subj-start':l[8],'subj-stop':l[9],'id':l[2],'Evalue':l[10],'bitScore':l[11],'len':l[3],'sens':stranded }]

    
#out.write("contig \t start \t stop \t subject \t sens \t id \t Evalue \t len \n") ## en tete
for listHit in bestHits: # parcourt de tous les scaf ayant un ou des eve
    i = 0 # indice pour la boucle while
    while i < len(bestHits[listHit]): # chaque eve du scaf
        out.write(bestHits[listHit][i]["contig"])
        out.write("\t")
        out.write(bestHits[listHit][i]["start"])
        out.write("\t")
        out.write(bestHits[listHit][i]["stop"])
        out.write("\t")
        out.write(bestHits[listHit][i]["subj"])
        out.write("\t")
        out.write(bestHits[listHit][i]["sens"])
        out.write("\t")
        out.write(bestHits[listHit][i]["id"])
        out.write("\t")
        out.write(bestHits[listHit][i]["Evalue"])
        out.write("\t")
        out.write(bestHits[listHit][i]["len"])
        out.write("\n")
        i+=1

# pprint.pprint(bestHits)
blast.close()
os.system("rm "+inFile ) ##
out.close()
