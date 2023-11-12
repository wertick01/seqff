from setuptools import setup, find_packages


def readme():
  with open('README.md', 'r') as f:
    return f.read()


setup(
  name='seqff',
  version='2.0.0',
  author='wertick01',
  author_email='cool.pumba01@yandex.ru',
  description='Fetal fraction calculating module/',
  long_description=readme(),
  long_description_content_type='text/markdown',
  url='https://github.com/wertick01',
  packages=find_packages(),
  install_requires=['matplotlib==3.5.1, numpy==1.26.1, pandas==1.4.2, pysam==0.21.0, rpy2==3.5.14, scikit_learn==1.0.2, scipy==1.7.3'],
  classifiers=[
    'Programming Language :: Python :: 3.10',
    # 'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent'
  ],
  keywords='files speedfiles ',
  project_urls={
    'GitHub': 'https://github.com/wertick01'
  },
  python_requires='>=3.6'
)