import collectionnames
import copy
import lciotools as lc
import math_utils as um
import mcana as mc

#get vertex number of a given rp in an event
#0 for ip, -ve for not in a vertex
def vertexNumber(rp, decay_chain_rptracks):
    #find the copy of the RP in the decay chain that points at our rp - if we can't then its isolated
    dcrp = [dcrp for dcrp in decay_chain_rptracks if dcrp.getParticles()[0].id() == rp.id()]
    #if itwasn't found its an isolated track
    if not dcrp: return -1
    #we now find the vertex number - trace back until we find a vertex with the primary flag set
    vertex = dcrp[0].getStartVertex()
    vertex_number = 0
    while not vertex.isPrimary():
        vertex_number += 1
        vertex = vertex.getAssociatedParticle().getStartVertex()
    return vertex_number
    
def getCorrectParent(mcp):
    pars = mcp.getParents()
    if len(pars) > 1:
        #pick the parent which is furthest from the IP if same return the first
        (flen,furthest) = (-1,0)
        for par in pars:
            if um.threeDRadius(par.getVertex()) > flen:
                (flen,furthest) = (um.threeDRadius(par.getVertex()),par)
        return furthest
    else:
        #An expection here means we made it to the root of the event
        #ie the ip detection below failed
        return pars[0]

def flavourOfParentVertex(event,rp):
    #note cannot just check direct parent due to resonances etc.
    track = rp.getTracks()[0]
    mcps = event.getRelatedTo(track,[collectionnames.trackmcpcollection,])
    mcp = mcps[0]
    mcp = getCorrectParent(mcp)
    #need to find the first parent that has a length unless we reach the IP (-1 quark type)
    while mc.quarkType(mcp.getPDG()) != -1 and mcp.getVertex() == mcp.getEndpoint():
        mcp = getCorrectParent(mcp)
    return mc.quarkType(mcp.getPDG())

def mcVertexNumber(event,rp):
    track = rp.getTracks()[0]
    mcps = event.getRelatedTo(track,[collectionnames.trackmcpcollection,])
    mcp = mcps[0]
    v = mcp.getVertex()
    #Go back through the tree counting all none resonant decays
    v_num = 0
    while mc.quarkType(mcp.getPDG()) != -1:
        #print mcp.getPDG()
        mcp = getCorrectParent(mcp)
        if um.threeDRadius(um.sub(mcp.getVertex(),mcp.getEndpoint())) > 0.0000000001 : v_num += 1
    #print flavourOfParentVertex(event,rp) , v_num, v
    return v_num
    
def numberOfVertices(decay_chain):
    #Have to compare by id
    return len(set([rp.getStartVertex().id() for rp in decay_chain.getParticles()]))
    
class vertexTable:
    def __init__(self,input_jet_rp_collection,decay_chain_collection, table_type="Parent_Flavour"):
        self.input_jet_rp_collection = input_jet_rp_collection
        self.decay_chain_collection = decay_chain_collection
        self.table_type = table_type
        self.results_table = {1:{"UDS": [0,0], "B":[0,0], "D":[0,0]},
                               2:{"UDS": [0,0,0], "B":[0,0,0], "D":[0,0,0]},
                               3:{"UDS": [0,0,0,0], "B":[0,0,0,0], "D":[0,0,0,0]}}
        self.num_entries = [0,0,0]

    def addJet(self,event,jet):
        dc = event.getCorrespondingObject(jet,self.decay_chain_collection) 
        num_verts= numberOfVertices(dc)
        if num_verts == 0: return
        if num_verts > 3: num_verts = 3
        self.num_entries[num_verts-1] = 1+self.num_entries[num_verts-1]
        if self.table_type == "Parent_Flavour":
            def flav(flav):
                if   flav == 4:return "D"
                elif flav == 5:return "B"
                else:          return "UDS"
            parent_flavours = [flav(flavourOfParentVertex(event,rp)) for rp in jet.getParticles()]
        elif self.table_type == "MCVNum":
            def vName(v_num):
                if   v_num == 0:return "UDS"
                elif v_num == 1:return "B"
                else:          return "D"
            parent_flavours = [vName(mcVertexNumber(event,rp)) for rp in jet.getParticles()]
        vertex_numbers = [vertexNumber(rp,dc.getParticles()) for rp in jet.getParticles()]
        for (p_flav,v_num) in zip(parent_flavours,vertex_numbers):
            if v_num > num_verts-1: v_num = num_verts-1
            #re-order so that iso at the end
            if v_num < 0: v_num = num_verts
            #print num_verts,p_flav,v_num
            self.results_table[num_verts][p_flav][v_num] += 1
                
    def printNormalisedTable(self):
        print sum(self.num_entries), "jets"
        for ((verts,sub_table),num_entries) in zip(self.results_table.iteritems(),self.num_entries):
            if num_entries:
                print verts, "vertices", '%.2f' % (float(num_entries)/float(sum(self.num_entries))*100.0)
                for (flav,row) in sub_table.iteritems():
                    if sum(row) > 0: print flav,"\t",['%.2f' % (float(entry)/float(sum(row))*100) for entry in row]                
    
    def normalisedTable(self):
        return [[[int(float(entry)/float(sum(row))*100) for entry in row] for (flav,row) in sub_table.iteritems()] for (verts,sub_table) in self.results_table.iteritems()]
          
    def printDiffTo(self, other):
        return [[[m-d for (m,d) in zip(rowm,rowd)]for (rowm,rowd) in zip(subtm,subtd)]for (subtm,subtd) in zip(self.normalisedTable(),other.normalisedTable())]

def vertexTables(event_gen, input_jet_rp_collection, decay_chain_collection, table_type = "Parent_Flavour"):
    tables = {}
    for name in ["ALL", "B", "D", "UDS"]:
        tables[name] = vertexTable(input_jet_rp_collection, decay_chain_collection, table_type)
    for e in event_gen:
        for jet in e[input_jet_rp_collection]:
            tables["ALL"].addJet(e,jet)
            true_flavour = e.getCorrespondingObject(jet,"TrueJetFlavour")[0]
            if true_flavour == 5:
                tables["B"].addJet(e,jet)
            elif true_flavour == 4:
                tables["D"].addJet(e,jet)
            else:
                tables["UDS"].addJet(e,jet)
    return tables
