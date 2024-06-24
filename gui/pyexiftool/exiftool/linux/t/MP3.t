# Before "make install", this script should be runnable with "make test".
# After "make install" it should work as "perl t/MP3.t".

BEGIN {
    $| = 1; print "1..2\n"; $Image::ExifTool::configFile = '';
    require './t/TestLib.pm'; t::TestLib->import();
}
END {print "not ok 1\n" unless $loaded;}

# test 1: Load the module(s)
use Image::ExifTool 'ImageInfo';
use Image::ExifTool::ID3;
use Image::ExifTool::MPEG;
$loaded = 1;
print "ok 1\n";

my $testname = 'MP3';
my $testnum = 1;

# test 2: Extract information from test image
{
    ++$testnum;
    my $exifTool = Image::ExifTool->new;
    my $info = $exifTool->ImageInfo('t/images/MP3.mp3');
    notOK() unless check($exifTool, $info, $testname, $testnum);
    print "ok $testnum\n";
}

done(); # end
