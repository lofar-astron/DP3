import os
import time

class Parset(object):

    def __init__(self, defaults = None):
        self.parameters = defaults or {}

    def readFromFile(self, fileName):
        lastline = ''
        for line in open(fileName, 'r').readlines():
            lastline = lastline + line.split('#')[0]
            lastline = lastline.rstrip()
            if len(lastline) > 0 and lastline[-1] == '\\':
                lastline = lastline[:-1]
            elif '=' in lastline:
                key, value = lastline.split('=')
                self.parameters[key.strip()] = value.strip()
                lastline = ''

    def writeToFile(self, fileName):
        outf = open(fileName, 'w')
        for key, value in sorted(self.parameters.iteritems()):
            outf.write(key + ' = ' + str(value) + '\n')
        outf.close()

    def getString(self, key):
        return self.parameters[key]

    def getInt32(self, key):
        return int(self.parameters[key])

    def getFloat(self, key):
        return float(self.parameters[key])

    def getStringVector(self, key):
        line = self.parameters[key]
        line.strip('[').strip(']')
        # this doesn't support \" in the string
        return line.split(',')

    def getInt32Vector(self, key):
        line = self.parameters[key]
        line.strip('[').strip(']')
        return [int(lp) for lp in line.split(',')]

    def getFloatVector(self, key):
        line = self.parameters[key]
        line.strip('[').strip(']')
        return [float(lp) for lp in line.split(',')]

    def __contains__(self, key):
        return key in self.parameters
        
    def __setitem__(self, key, value):
        self.parameters[key] = value

    def __getitem__(self, key):
        return self.parameters[key]



