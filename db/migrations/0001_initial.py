# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'SmallMolecule'
        db.create_table(u'db_smallmolecule', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sm_id', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8)),
            ('name', self.gf('django.db.models.fields.TextField')()),
            ('smiles', self.gf('django.db.models.fields.TextField')()),
            ('inchi', self.gf('django.db.models.fields.TextField')()),
            ('inchi_key', self.gf('django.db.models.fields.TextField')()),
            ('IUPAC_Name', self.gf('django.db.models.fields.TextField')()),
            ('alternative_names', self.gf('django.db.models.fields.TextField')()),
            ('Metabolic_Network_ids', self.gf('django.db.models.fields.TextField')()),
            ('pubchem_cid', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('chebi_id', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('kegg_id', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('hmdb_id', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('drugbank_id', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal(u'db', ['SmallMolecule'])

        # Adding model 'Gene'
        db.create_table(u'db_gene', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('gene_id', self.gf('django.db.models.fields.CharField')(unique=True, max_length=8)),
            ('name', self.gf('django.db.models.fields.TextField')()),
            ('alternate_name', self.gf('django.db.models.fields.TextField')()),
            ('entrez_id', self.gf('django.db.models.fields.TextField')()),
            ('entrez_symbol', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['Gene'])


    def backwards(self, orm):
        # Deleting model 'SmallMolecule'
        db.delete_table(u'db_smallmolecule')

        # Deleting model 'Gene'
        db.delete_table(u'db_gene')


    models = {
        u'db.gene': {
            'Meta': {'object_name': 'Gene'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'entrez_id': ('django.db.models.fields.TextField', [], {}),
            'entrez_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene_id': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '8'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {})
        },
        u'db.smallmolecule': {
            'IUPAC_Name': ('django.db.models.fields.TextField', [], {}),
            'Meta': {'object_name': 'SmallMolecule'},
            'Metabolic_Network_ids': ('django.db.models.fields.TextField', [], {}),
            'alternative_names': ('django.db.models.fields.TextField', [], {}),
            'chebi_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'drugbank_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'hmdb_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'inchi': ('django.db.models.fields.TextField', [], {}),
            'inchi_key': ('django.db.models.fields.TextField', [], {}),
            'kegg_id': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'pubchem_cid': ('django.db.models.fields.CharField', [], {'max_length': '15'}),
            'sm_id': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '8'}),
            'smiles': ('django.db.models.fields.TextField', [], {})
        }
    }

    complete_apps = ['db']