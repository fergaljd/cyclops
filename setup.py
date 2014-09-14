#!/usr/bin/env python

from distutils.core import setup


setup(name='PepLibGen',
      version='1.2',
      url='',
      author_email='fergaljd@gmail.ie',
      packages=['PepLibGen', 'PepLibGen.StructGen', 'PepLibGen.Analysis'],
      scripts = [
              'CycloPs/CycloPs.py'
          ]
      )