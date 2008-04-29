import lcio
from itertools import chain
import os
import pdgnames

def vecToList(vec):
    return_list = [] 
    for i in range(vec.size()):
        return_list.append(vec[i])
    return return_list

def listFrom(collection):
    return [collection.getElementAt(i) for i in xrange(collection.getNumberOfElements())]
   
def genFrom(collection):
    return (collection.getElementAt(i) for i in xrange(collection.getNumberOfElements()))
    
def getReader(file):
    fac=lcio.LCFactory.getInstance()
    rdr=fac.createLCReader()
    rdr.open(file)
    return rdr

class Event:
    def __init__(self,lcio_evt): 
        self.lcio_evt = lcio_evt
        self.relatedTo = {}
        self.collections= {}
        
    def __getitem__(self,collection_name):
        try:
            return (o for o in self.collections[collection_name])
        except KeyError:
            self.collections[collection_name] = listFrom(self.lcio_evt.getCollection(collection_name))
            return (o for o in self.collections[collection_name])
    
    def collectionNames(self):
        return vecToList(self.lcio_evt.getCollectionNames())
    
    def getCollection(self,collection_name):
        return self.lcio_evt.getCollection(collection_name)
        
    def getObjectByID(self,lc_object_id):
        try:
            return self.objects[lc_object_id]
        except AttributeError:
            collections = [self.lcio_evt.getCollection(name) for name in self.collectionNames()]
            objects = [listFrom(collection) for collection in collections]
            objects = chain(*objects)
            self.objects = {}
            for obj in objects: self.objects[obj.id()] = obj
            return self.objects[lc_object_id]
        
    def getRelatedTo(self,lc_object,collection_names=None):
        if not collection_names:
           collection_names = [name for name in self.collectionNames()
                               if self.getCollection(name).getTypeName() == "LCRelation"]
        try:
            return self.relatedTo[(lc_object,tuple(collection_names))]
        except KeyError:
            relations = [self[collection]
                         for collection in collection_names]
            #flatten the list
            relations = list(chain(*relations))
            #We need to get relations with our object in both the from and to:
            related_from = [relation.getFrom()
                            for relation in relations
                            if relation.getTo().id() == lc_object.id()]
            related_to = [relation.getTo()
                          for relation in relations 
                          if relation.getFrom().id() == lc_object.id()]
            related = related_from + related_to
            #We have LCObjects - would be better if they were the sub-class
            related = [self.getObjectByID(obj.id())
                       for obj in related]
            self.relatedTo[(lc_object,tuple(collection_names))] = related
            return related

    def getObjFromID(self,idn, collection):
        objs = list(self[collection])
        ids = [obj.id() for obj in objs]
        return objs[ids.index(idn)]
                
    def getCorrespondingObject(self,original_object,corresponding_collection):
        corresponding_collection = list(self[corresponding_collection])
        #find which collection the original_object is in
        collections = [genFrom(self.getCollection(name))
                       for name in self.collectionNames()
                       if self.getCollection(name).getNumberOfElements() == len(corresponding_collection)]
        is_in_collection = [original_object.id() in [obj.id() for obj in collection] for collection in collections]
        collections = [self[name]
                       for name in self.collectionNames()
                       if self.getCollection(name).getNumberOfElements() == len(corresponding_collection)]
        #If this fails then no collection that contained the object    
        containing_collection = list(collections[is_in_collection.index(True)])
        #we should check that the two are the same length - otherwise its likely that what you asked for is not a corresponding var
        return corresponding_collection[[obj.id() for obj in containing_collection].index(original_object.id())]
        
def genEventsFromFile(f):
    madeReader = False
    while True:
        try:
            lcio_evt = f.readNextEvent()
        except AttributeError:
            f = getReader(f)
            madeReader = True
            lcio_evt = f.readNextEvent()
        if lcio_evt == None:
            if madeReader: f.close()
            return
        evt = Event(lcio_evt)
        yield evt

def genEventsFromDir(d):
    return chain(*(genEventsFromFile(os.path.join(d,f))
                 for f in os.listdir(d) 
                 if f.endswith(".slcio")))

