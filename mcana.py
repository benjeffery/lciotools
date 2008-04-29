import lciotools as lc
import math_utils as mu
import exceptions
import subprocess as sp
from copy import copy

class BadJetFinding(exceptions.Exception):
    def __init__(self,msg):
		return
    def __str__(self):
		print "",msg

def quarkType(pdg):
    if pdg in (91,92,93): return -1
    pdg = str(abs(pdg))
    if len(pdg) == 1: return int(pdg)
    if len(pdg) == 2: return 0
    if len(pdg) == 4: return int(pdg[:1])
    return int(pdg[-3])
    
def getTrueJetFlavour(event, jet):
    return lc.pdgToName(event.getCorrespondingObject(jet,"TrueJetFlavour")[0])
    
def getMCHeavyPartonOf(event,jet):
    #For each particle in the jet get its MCP and trace that back to its initial parton
    jet_trks = [p.getTracks() for p in jet.getParticles()]
    jet_trks = sum(jet_trks,())
    mcps = [event.getRelatedTo(trk) for trk in jet_trks]
    mcps = sum(mcps,[])
    
    def findParton(mcp):
        last = mcp
        while mcp.getParents() and mcp.getParents()[0].getPDG() != 92:
            mcp = mcp.getParents()[0]
            last = mcp
        return last
    
    partons = [findParton(mcp) for mcp in mcps]
    partons = [p.id() for p in partons
               if quarkType(p.getPDG()) > 3]
    #return partons
    unique_set = set(partons)
    if len(unique_set) > 1: raise BadJetFinding("More than one parton associated with jet")
    #convert id's back to mcps
    try:
        return event.getObjFromID(unique_set.pop(),"MCParticle")
    except KeyError:
        return None

def getHeavyHadronDecayLength(mcp):
    #True if flavour change
    start = mcp.getVertex()
    same_flavour = [quarkType(d.getPDG()) == quarkType(mcp.getPDG()) for d in mcp.getDaughters()]
    while any(same_flavour):
        mcp = mcp.getDaughters()[same_flavour.index(True)]
        same_flavour = [quarkType(d.getPDG()) == quarkType(mcp.getPDG()) for d in mcp.getDaughters()]
    end = mcp.getEndpoint()
    return mu.threeDRadius(mu.sub(start,end))
    
    
def mCVertices(mcps):
    #group mc particles by production point, assume same point if they have a distance less than 0.01micron!
    vertices = {}
    for mcp in mcps:
        #check that there is not a key that is below the threshold
        nearkeys = (key
                    for key in vertices.iterkeys()
                    if mu.threeDRadius(mu.sub(mcp.getVertex(),key)) < 0.0000001)
        try:
            vertices[nearkeys.next()].append(mcp)
        except StopIteration:
            try:
                vertices[mcp.getVertex()].append(mcp)
            except KeyError:
                vertices[mcp.getVertex()] = [mcp,]
    return vertices

def fromIP(mcp):
    if mcp.getPDG() == 92: return True
    while mcp.getParents():
        if mcp.getParents()[0].getPDG() == 92: 
            return True
        mcp = mcp.getParents()[0]
    return False

def fromHeavy(mcp):
    if quarkType(mcp.getPDG()) > 3: return True
    if mcp.getPDG() == 92: return True
    while mcp.getParents() and mcp.getPDG() != 92:
        if quarkType(mcp.getParents()[0].getPDG()) > 3: 
            return True
        mcp = mcp.getParents()[0]
    return False


def printMCTree(mcps,filename="temp.dot"):
    #start the graph
    f=open(filename, 'w')
    f.write('digraph G {\nranksep="equally";\noverlap="false";\nrankdir="LR";\ncompound=true;\n')
    #make a node for each mcp 
    colors = {-1:"white",0:"white",1:"grey",2:"yellow",3:"green",4:"red",5:"blue",6:"pink"}
    for mcp in mcps:
        name = str(mcp.id())
        label = lc.pdgToName(mcp.getPDG())
        colour= colors[quarkType(mcp.getPDG())]
        f.write('"'+name+'" [ label="'+label+' '+str(len(mcp.getParents()))+'",style="filled",color="'+colour+'" ];\n')
    #idn = 0
    #mcp_ids = [mcp.id() for mcp in mcps]
    #for mcp in mcps:
    #    daughters = [daughter for daughter in mcp.getDaughters() if daughter.id() in mcp_ids]
    #    if daughters:
    #        f.write('subgraph cluster'+str(idn)+'{\n')
    #        for daughter in daughters:
    #            f.write(str(daughter.id())+';\n')
    #        f.write('}\n')
    #        f.write(str(mcp.id())+'->'+str(daughter.id())+' [ lhead=cluster'+str(idn)+',label = #"'+str(mu.threeDRadius(sub(mcp.getVertex(),mcp.getEndpoint())))[:4]+'" ];')
    #        idn += 1

    #make a node in each cluster
    #idn = 0
    #for vert in mCVertices(mcps).itervalues():
    #    #we need to remove the ones that leave this vertex
    #    #vert = [mcp for mcp in vert if mu.threeDRadius(sub(mcp.getVertex(),mcp.getEndpoint())) < 0.00000001]
    #    if len(vert) > 1:
    #        f.write('subgraph cluster'+str(idn)+'{\n')
    #        for mcp in vert:
    #            f.write(str(mcp.id())+';\n')
    #        f.write('}\n')
    #        idn += 1
    #make a link for each decay (whose product is in the list!)
    mcp_ids = [mcp.id() for mcp in mcps ]#if fromIP(mcp)]
    for mcp in mcps:
        for parent in mcp.getParents():
     #       if fromIP(parent):
                f.write(str(parent.id())+'->'+str(mcp.id())+' [ label = "'+str(mu.threeDRadius(mu.sub(parent.getVertex(),parent.getEndpoint())))[:4]+'" ];')
    #end graph
    f.write('}\n')
    f.close()
    f=open(filename+".png", 'w')
    sp.call(['dot', '-Tpng', filename],stdout = f)
    f.close()
    sp.call(['eog', filename+".png"])


