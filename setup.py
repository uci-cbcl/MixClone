from distutils.core import setup

setup(
      name = 'MixClone',
      version = '1.0',
      description = 'MixClone: a mixture model for inferring tumor subclonal populations',
      author = 'Yi Li, Andrew Roth',
      author_email = 'yil8@uci.edu, andrewjlroth@gmail.com',
      url = 'https://github.com/uci-cbcl/MixClone',
      license = 'GNU GPL v2',
      packages = [
                  'mixclone',
                  'mixclone.preprocess',
                  'mixclone.model',
                  'mixclone.postprocess'
      ],
      scripts = ['MixClone.py']
)
