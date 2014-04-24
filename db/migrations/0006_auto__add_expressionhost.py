# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'ExpressionHost'
        db.create_table(u'db_expressionhost', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('host_id', self.gf('django.db.models.fields.CharField')(default='UWPRHWYT', unique=True, max_length=8)),
            ('name', self.gf('django.db.models.fields.TextField')()),
            ('taxonomic_id', self.gf('django.db.models.fields.TextField')()),
            ('taxonomic_branch', self.gf('django.db.models.fields.TextField')()),
            ('genotype', self.gf('django.db.models.fields.CharField')(max_length=32)),
            ('genome', self.gf('django.db.models.fields.TextField')()),
            ('expression_data', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('metabolic_model', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('reference', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['ExpressionHost'])


    def backwards(self, orm):
        # Deleting model 'ExpressionHost'
        db.delete_table(u'db_expressionhost')


    models = {
        u'db.expressionhost': {
            'Meta': {'object_name': 'ExpressionHost'},
            'expression_data': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'genome': ('django.db.models.fields.TextField', [], {}),
            'genotype': ('django.db.models.fields.CharField', [], {'max_length': '32'}),
            'host_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRHWZB'", 'unique': 'True', 'max_length': '8'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'metabolic_model': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'reference': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'taxonomic_branch': ('django.db.models.fields.TextField', [], {}),
            'taxonomic_id': ('django.db.models.fields.TextField', [], {})
        },
        u'db.gene': {
            'Meta': {'object_name': 'Gene'},
            'alternate_name': ('django.db.models.fields.TextField', [], {}),
            'entrez_id': ('django.db.models.fields.TextField', [], {}),
            'entrez_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRHWZF'", 'unique': 'True', 'max_length': '8'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {})
        },
        u'db.protein': {
            'Meta': {'object_name': 'Protein'},
            'amino_acid_sequence': ('django.db.models.fields.TextField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'protein_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRHWZA'", 'unique': 'True', 'max_length': '8'}),
            'uniprot_id': ('django.db.models.fields.CharField', [], {'max_length': '15'})
        },
        u'db.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'alternate_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'chemical_equation': ('django.db.models.fields.TextField', [], {}),
            'comment': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'confidence_score': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'directionality': ('django.db.models.fields.TextField', [], {}),
            'ec_number': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kegg_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'reaction_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRHWZD'", 'unique': 'True', 'max_length': '8'}),
            'reference': ('django.db.models.fields.TextField', [], {'blank': 'True'})
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
            'sm_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRHWZE'", 'unique': 'True', 'max_length': '8'}),
            'smiles': ('django.db.models.fields.TextField', [], {})
        },
        u'db.substance': {
            'Meta': {'object_name': 'Substance'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        }
    }

    complete_apps = ['db']