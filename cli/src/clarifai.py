#!/usr/bin/env python

import os
import sys

from clarifai.client import ClarifaiApi


def tag_images_in_directory(path, api):
  images = []
  path = path.rstrip(os.sep)
  for fname in os.listdir(path):
    images.append((open(os.path.join(path, fname), 'rb'), fname))
  return api.tag_images(images)


def main(argv):
  if len(argv) > 1:
    imageurl = argv[1]
  else:
    imageurl = 'http://clarifai-img.s3.amazonaws.com/test/toddler-flowers.jpeg'

  api = ClarifaiApi()

  if imageurl.startswith('http'):
    response = api.tag_image_urls(imageurl)
  elif os.path.isdir(imageurl):
    response = tag_images_in_directory(imageurl, api)
  elif os.path.isfile(imageurl):
    with open(imageurl,'rb') as image_file:
      response = api.tag_images(image_file)
  else:
    raise Exception("Must input url, directory path, or file path")

  print(response)


if __name__ == '__main__':
  main(sys.argv)
