# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Protein'
        db.create_table(u'db_protein', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('protein_id', self.gf('django.db.models.fields.CharField')(default='UWPRGWV2', unique=True, max_length=8)),
            ('uniprot_id', self.gf('django.db.models.fields.CharField')(max_length=15)),
            ('amino_acid_sequence', self.gf('django.db.models.fields.TextField')()),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Gene'])),
        ))
        db.send_create_signal(u'db', ['Protein'])

        # Adding field 'Gene.sequence'
        db.add_column(u'db_gene', 'sequence',
                      self.gf('django.db.models.fields.TextField')(default=''),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting model 'Protein'
        db.delete_table(u'db_protein')

        # Deleting field 'Gene.sequence'
        db.delete_column(u'db_gene', 'sequence')


    models = {
        u'db.gene': {
            'Meta': {'object_name': 'Gene'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'entrez_id': ('django.db.models.fields.TextField', [], {}),
            'entrez_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWV9'", 'unique': 'True', 'max_length': '8'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {})
        },
        u'db.protein': {
            'Meta': {'object_name': 'Protein'},
            'amino_acid_sequence': ('django.db.models.fields.TextField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWV8'", 'unique': 'True', 'max_length': '8'}),
            'uniprot_id': ('django.db.models.fields.CharField', [], {'max_length': '15'})
        },
        u'db.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'chemical_equation': ('django.db.models.fields.TextField', [], {}),
            'comment': ('django.db.models.fields.TextField', [], {}),
            'confidence_score': ('django.db.models.fields.TextField', [], {}),
            'directionality': ('django.db.models.fields.TextField', [], {}),
            'ec_number': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kegg_id': ('django.db.models.fields.TextField', [], {}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'reaction_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWV7'", 'unique': 'True', 'max_length': '8'}),
            'reference': ('django.db.models.fields.TextField', [], {})
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
            'sm_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRGWVA'", 'unique': 'True', 'max_length': '8'}),
            'smiles': ('django.db.models.fields.TextField', [], {})
        },
        u'db.substance': {
            'Meta': {'object_name': 'Substance'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        }
    }

    complete_apps = ['db']