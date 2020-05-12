# -*- coding: utf-8 -*-

from __future__ import unicode_literals

import unittest
import exiftool
import warnings
import os

class TestExifTool(unittest.TestCase):
    def setUp(self):
        self.et = exiftool.ExifTool()
    def tearDown(self):
        if hasattr(self, "et"):
            self.et.terminate()
        if hasattr(self, "process"):
            if self.process.poll() is None:
                self.process.terminate()
    def test_termination_cm(self):
        # Test correct subprocess start and termination when using
        # self.et as a context manager
        self.assertFalse(self.et.running)
        self.assertRaises(ValueError, self.et.execute)
        with self.et:
            self.assertTrue(self.et.running)
            with warnings.catch_warnings(record=True) as w:
                self.et.start()
                self.assertEquals(len(w), 1)
                self.assertTrue(issubclass(w[0].category, UserWarning))
            self.process = self.et._process
            self.assertEqual(self.process.poll(), None)
        self.assertFalse(self.et.running)
        self.assertNotEqual(self.process.poll(), None)
    def test_termination_explicit(self):
        # Test correct subprocess start and termination when
        # explicitly using start() and terminate()
        self.et.start()
        self.process = self.et._process
        self.assertEqual(self.process.poll(), None)
        self.et.terminate()
        self.assertNotEqual(self.process.poll(), None)
    def test_termination_implicit(self):
        # Test implicit process termination on garbage collection
        self.et.start()
        self.process = self.et._process
        del self.et
        self.assertNotEqual(self.process.poll(), None)
    def test_get_metadata(self):
        expected_data = [{"SourceFile": "rose.jpg",
                          "File:FileType": "JPEG",
                          "File:ImageWidth": 70,
                          "File:ImageHeight": 46,
                          "XMP:Subject": "Röschen",
                          "Composite:ImageSize": "70x46"},
                         {"SourceFile": "skyblue.png",
                          "File:FileType": "PNG",
                          "PNG:ImageWidth": 64,
                          "PNG:ImageHeight": 64,
                          "Composite:ImageSize": "64x64"}]
        script_path = os.path.dirname(__file__)
        source_files = []
        for d in expected_data:
            d["SourceFile"] = f = os.path.join(script_path, d["SourceFile"])
            self.assertTrue(os.path.exists(f))
            source_files.append(f)
        with self.et:
            actual_data = self.et.get_metadata_batch(source_files)
            tags0 = self.et.get_tags(["XMP:Subject"], source_files[0])
            tag0 = self.et.get_tag("XMP:Subject", source_files[0])
        for expected, actual in zip(expected_data, actual_data):
            et_version = actual["ExifTool:ExifToolVersion"]
            self.assertTrue(isinstance(et_version, float))
            if isinstance(et_version, float):    # avoid exception in Py3k
                self.assertTrue(
                    et_version >= 8.40,
                    "you should at least use ExifTool version 8.40")
            actual["SourceFile"] = os.path.normpath(actual["SourceFile"])
            for k, v in expected.items():
                self.assertEqual(actual[k], v)
        tags0["SourceFile"] = os.path.normpath(tags0["SourceFile"])
        self.assertEqual(tags0, dict((k, expected_data[0][k])
                                     for k in ["SourceFile", "XMP:Subject"]))
        self.assertEqual(tag0, "Röschen")

if __name__ == '__main__':
    unittest.main()
