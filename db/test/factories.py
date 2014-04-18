import factory
import db.models

class SmallMoleculeFactory(factory.django.DjangoModelFactory):
    FACTORY_FOR=db.models.SmallMolecule

    sm_id = factory.Sequence(lambda n: str(n))
    name = 'name_%s' % str(factory.Sequence(lambda n: str(n)))
    smiles = ('CC[C@H](CO)NC1=NC2=C(C(=N1)NCC3=CC=CC=C3)N=CN2C(C)C%s' 
                % str(factory.Sequence(lambda n: str(n))) )
