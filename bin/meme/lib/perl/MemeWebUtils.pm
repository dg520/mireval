# File: MemeWebUtils.pm
# Project: Website CGI
# Description: MemeWebUtils.pm made from MemeWebUtils.pm.in by make. Helper functions for CGI pages.

package MemeWebUtils;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(is_safe_name get_safe_name get_safe_file_name getnum 
    is_numeric get_alph valid_address valid_meme_version add_status_msg 
    update_status loggable_date write_invocation_log eps_to_png);

use Cwd;
use Fcntl qw(O_APPEND O_CREAT O_WRONLY O_TRUNC SEEK_SET);
use File::Basename qw(fileparse);
use File::Copy qw(copy);
use File::Spec::Functions qw(catfile splitdir tmpdir);
use File::Temp qw(tempfile);
use HTTP::Request::Common qw(GET);
use XML::Simple;
use HTML::PullParser;
use HTML::Template;
use Scalar::Util qw(openhandle);
use Sys::Hostname;

use lib qw(/apps/lab/miket/miREval/2.0/mireval/bin/meme/lib/perl);
use CatList qw(load_categories load_entry);
use ExecUtils qw(invoke);

# Setup logging
my $logger = undef;
eval {
  require Log::Log4perl;
  Log::Log4perl->import();
};
unless ($@) {
  Log::Log4perl::init('/apps/lab/miket/miREval/2.0/mireval/bin/meme/etc/logging.conf');
  $logger = Log::Log4perl->get_logger('meme.cgi.utils');
}

my $MOTIF_DB_INDEX = '/apps/lab/miket/miREval/2.0/mireval/bin/meme/etc/motif_db.index';
my $template_dir = '/apps/lab/miket/miREval/2.0/mireval/bin/meme/web/cgi-bin';
my $service_invocation_log_dir = '/apps/lab/miket/miREval/2.0/mireval/bin/meme/LOGS';
my $tmpdir = '';
# use the perl default if none is supplied or the replace fails
$tmpdir = &tmpdir() if ($tmpdir eq '' || $tmpdir =~ m/^\@TMP[_]DIR\@$/);
my @gs_version_nums = ();
my $GHOSTSCRIPT = '/usr/bin/gs';
my $CONVERT = '/usr/bin/convert';
##############################################################################
#          Functions
##############################################################################

sub log_and_die {
  if ($logger) {
    $Log::Log4perl::caller_depth++;
    $logger->loglog_and_die(@_);
  } else {
    die(@_);
  }
}

my $SAFE_NAME_CHARACTERS = 'a-zA-Z0-9:_\.\-';

# Checks that a file name has only whitelisted characters in it and does
# not have a leading dash
sub is_safe_name {
  $logger->trace('call is_safe_name') if $logger;
  my ($name) = @_;
  if ($name =~ /^[$SAFE_NAME_CHARACTERS]*$/ && $name !~ /^-/) {
    return 1;
  }
  return 0;
}

# Removes characters from a string that could cause problems in the shell
# Note that it is possible that the name is reduced to nothing so you can
# specify an alternate name and the length that the safe name must
# exceed to avoid using the alternate name.
sub get_safe_name {
  $logger->trace('call get_safe_name') if $logger;
  my ($name, $alt, $len) = @_;
  return $alt unless $name;
  $len = 0 unless defined($len);
  $name =~ tr/ /_/;
  $name =~ s/[^$SAFE_NAME_CHARACTERS]//g;
  #remove leading dashes
  $name =~ s/^-+//g;
  #this is primarily for satisfying taint mode
  if ($name =~ /^([$SAFE_NAME_CHARACTERS]*)$/) { 
    $name = $1;
  } else { #should never be executed
    log_and_die("Name contains invalid characters");
  }
  if (length($name) > $len) {
    return $name;
  } else {
    return $alt;
  }
}

# Removes any path and calls get_safe_name
sub get_safe_file_name {
  $logger->trace('call get_safe_file_name') if $logger;
  my ($file_name, $alt, $len) = @_;
  return $alt unless $file_name;
  $file_name = fileparse($file_name); # remove directory component
  return get_safe_name($file_name, $alt, $len);
}

#
# Use POSIX strtod to convert a variable to a number.
#
# Used in
# Utils
#
sub getnum {
  $logger->trace('call getnum') if $logger;

  use POSIX qw(strtod);
  my $str = shift;
  $str =~ s/^\s+//;
  $str =~ s/\s+$//;
  $! = 0;
  my($num, $unparsed) = strtod($str);
  if (($str eq '') || ($unparsed != 0) || $!) {
      return undef;
  } else {
      return $num;
  } 

} 

#
# Test whether a variable is a number.
#
# Used in
# Utils
#
sub is_numeric { 
  $logger->trace('call is_numeric') if $logger;
  defined &getnum($_[0]) 
} 

#
# Output large integers with commas
#
#
sub commify {
  $logger->trace('call commify') if $logger;
  my $input = shift;
  $input = reverse $input;
  $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return reverse $input;
} 

#
# get the alphabet of a string: DNA or PROTEIN
#
# Used in
# Utils
#
sub get_alph
{
  $logger->trace('call get_alph') if $logger;
  my ($letters) = @_;                # get arguments 
  my ($old, $new, $x);

  $old = length($letters);

  # check against allowed nucleotide letters
  $x = $letters;
  $x =~ tr/ABCDGHKMNRSTUVWY//cd; # delete all letters that aren't a nucleotide letter
  $new = length($x);
  if ($old == $new) { # if nothing was deleted then it has only dna letters
    return("DNA", "");
  } else {
    # check against allowed protein letters
    $x = $letters;
    $x =~ tr/ABCDEFGHIKLMNPQRSTUVWXYZ//cd; # delete all non-protein letters
    $new = length($x);
    if ($old == $new) { # if nothing was deleted then it has only protein letters
      return("PROTEIN", "");
    } else {
      # get the unknown letters
      $x = $letters;
      $x =~ tr/ABCDEFGHIKLMNPQRSTUVWXYZ//d; # delete all known letters
      return("UNRECOGNIZED", $x);
    }
  }
} # get_alph

#
# valid_address
#
# Checks email address syntax and host validity.
# Returns 0 if the address is invalid, 1 on success.
#
# Merged from Validation.txt
#
# Used in
# Utils
#
sub valid_address {
  $logger->trace('call valid_address') if $logger;
  my($addr) = @_;
  my($domain, $valid);
  return(0) unless ($addr =~ /^[^@]+@([-\w]+\.)+[A-Za-z]{2,4}$/);
  $domain = (split(/@/, $addr))[1];
  $valid = 0; open(DNS, "nslookup -q=mx $domain |") || return(-1);
  while (<DNS>) {
    #my $line = (/^$domain.*mail exchanger/);
    $valid = 1 if (/^$domain.*\s(mail exchanger|internet address)\s=/);
  }
  return($valid);
}

#
# valid_meme_version
#
# Checks a version string to ensure that it should be parsable.
# Return 0 if the version is unknown or unparsable, 1 on success.
#
sub valid_meme_version {
  $logger->trace('call valid_meme_version') if $logger;
  my ($ver_str) = @_;
  my @ver_strs = split(/\./, $ver_str);
  #check that there is a version number at all
  return 0 unless scalar(@ver_strs) > 0;
  #check that each part of a version number is actually a number
  my @ver_nums = ();
  foreach my $ver_part (@ver_strs) {
    my $ver_num = getnum($ver_part);
    return 0 unless defined $ver_num;
    push(@ver_nums, $ver_num);
  }
  # check the major version number
  if ($ver_nums[0] < 3) {
    return 0; # prior to version 3 is too old
  }
  # get the current version numbers
  my $cur_str = "4.9.1";
  my @cur_nums = split(/\./, $cur_str);
  # compare the current version to the specified
  while (@ver_nums && @cur_nums) {
    my $ver_num = shift(@ver_nums);
    my $cur_num = shift(@cur_nums);
    if ($ver_num > $cur_num) {
      # it's newer, hence unknown
      return 0;
    } elsif ($ver_num < $cur_num) {
      # older
      last;
    }
  }
  #seems to be valid
  return 1;
}


#
# add_status_msg
#
# Adds a message to the message list
# and returns the list. For use with
# update_status.
#
sub add_status_msg {
  $logger->trace('call add_status_msg') if $logger;
  my ($msg, $msg_list) = @_;
  push(@$msg_list, {msg => $msg});
  return $msg_list;
}

#
# update_status
#
# Creates or updates the specified status file to contain the
# file list only showing each of the files if they exist and
# the message list.
#
sub update_status {
  $logger->trace('call update_status') if $logger;
  my ($output_file, $program, $refresh, $files_list, $msg_list, $status) = @_;

  my @found_files = ();
  foreach my $entry (@$files_list) {
    my $file = $entry->{'file'};
    push(@found_files, $entry) if (defined($file) && -e $file && -s $file);
  }

  my $fh;
  sysopen($fh, $output_file, O_CREAT | O_WRONLY | O_TRUNC) or log_and_die("Failed to open \"$output_file\".");
  my $template = HTML::Template->new(filename => "$template_dir/job_status.tmpl");
  $template->param(program => $program, refresh => $refresh, files => \@found_files, msgs => $msg_list, status => $status);
  print $fh $template->output;
  close($fh) or log_and_die("Failed to close \"$output_file\".");
}

#
# loggable_date
#
# Creates a date string that is suitable for putting in the service invocation log.
#
sub loggable_date {
  $logger->trace('call loggable_date') if $logger;
  my $timestamp = shift;
  $timestamp = time() unless defined($timestamp);
  # old method 
  # return `date -u '+%d/%m/%y %H:%M:%S'`;
  # new method that doesn't require forking a new process
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime($timestamp);
  return sprintf('%02d/%02d/%02d %02d:%02d:%02d', $mday, $mon + 1, $year % 100, $hour, $min, $sec);
}

#
#
# write_invocation_log 
#
# Writes the date and time of a service's invocation to a log file
#
sub write_invocation_log {
  $logger->trace('call write_invocation_log') if $logger;
  my ($file, $start_time, $args) = @_;
  # the host
  my $host = hostname();
  # the current directory without path
  my $jobid = (splitdir(getcwd()))[-1];
  # the submission time if it is avaliable but use the start time as a default
  my $submit_time = $start_time;
  if (-e 'submit_time_file') {
    $submit_time = `cat submit_time_file`;
    unlink 'submit_time_file';
  }
  # the end time (now)
  my $end_time = loggable_date();
  # the email if it is avaliable but use our support address as a default
  my $email = 'meme@nbcr.net';
  if (-e 'address_file') {
    $email = `cat address_file`;
    unlink 'address_file';
  }

  my $logfile = catfile($service_invocation_log_dir, $file);
  my $logfh;
  sysopen($logfh, $logfile, O_CREAT | O_APPEND | O_WRONLY) or log_and_die("Unable to open invocation log for appending.");
  print $logfh "$host $jobid submit: $submit_time start: $start_time end: $end_time $args $email\n"; 
  close($logfh);

}

#
# eps_to_png
#
# Converts an EPS file into a png file
#
sub eps_to_png {
  $logger->trace('call eps_to_png') if $logger;
  my ($eps_filename, $png_filename, $messages) = @_;
  my @gs_args = ('-q', '-r100', '-dSAFER', '-dBATCH', '-dNOPAUSE', 
    '-dDOINTERPOLATE', '-sDEVICE=pngalpha', '-dBackgroundColor=16#FFFFFF', 
    '-dTextAlphaBits=4', '-dGraphicsAlphaBits=4', '-dEPSCrop', 
    '-sOutputFile='. $png_filename, $eps_filename);
  my @convert_args = ('eps:'.$eps_filename, 'png:'.$png_filename);
  my $status = -1;
  if (&_gs_ok()) {
    $status = invoke(
      PROG => $GHOSTSCRIPT,
      ARGS => \@gs_args,
      ALL_VAR => $messages
    );
  } elsif (&_im_ok()) {
    $status = invoke(
      PROG => $CONVERT,
      ARGS => \@convert_args,
      ALL_VAR => $messages
    );
    return system($CONVERT, @convert_args);
  } else {
    return 1;
  }
}

#
# Private Function
#
# Provide basic escaping of xml characters in a string.
# Note that does not detect existing escaped characters
# and so they would be re-escaped.
#
sub _fix_xml {
  $logger->trace('call fix_xml') if $logger;
  my $text = shift;
  $text =~ s/&/&amp;/g;
  $text =~ s/</&lt;/g;
  $text =~ s/>/&gt;/g;
  $text =~ s/"/&quot;/g;
  return $text;
}

#
# Private Function
#
# Gets the version of ghostscript installed
#
sub _gs_version {
  return @gs_version_nums if (@gs_version_nums);
  @gs_version_nums = (-1); # default to failure 
  if (-e $GHOSTSCRIPT && -x $GHOSTSCRIPT) {
    my @nums = map(int, split(/\./, `$GHOSTSCRIPT --version`));
    if (@nums) {
      @gs_version_nums = @nums;
    }  
  }
  return @gs_version_nums;
}

#
# Private Function
#
# Returns true if ghostscript is ok to use
#
sub _gs_ok {
  my @ver = _gs_version();
  return $ver[0] >= 8;
}

#
# Private Function
#
# Returns true if image magicks convert is ok to use
#
sub _im_ok {
  return -e $CONVERT && -x $CONVERT;
}

##############################################################################
#                      Object Methods
##############################################################################

#
# new
#
# Constructor.
#
# Used in
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub new {
  $logger->trace('call MemeWebUtils->new') if $logger;
  my $classname = shift;
  my $self = {};
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

#
# Private Method
#
# _init
#
# Sets the program to the passed value.
#
sub _init {
  $logger->trace('call MemeWebUtils->_init') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($program, $bin_dir) = @_;
  $program = "" unless defined $program;
  $self->{PROGRAM} = $program;
  $self->{BIN_DIR} = $bin_dir;
  $self->{NERRORS} = 0;
}

#
# get_program
#
# Gets the program name as passed to the constructor.
#
sub get_program {
  $logger->trace('call MemeWebUtils->get_program') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  return $self->{PROGRAM};
}

#
# has_errors
#
# Returns the number of errors that have been detected.
#
# Used in
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub has_errors {
  $logger->trace('call MemeWebUtils->has_errors') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  return $self->{NERRORS};
}

#
# Check to see whether the email address is valid.
#
# Used in
# dreme.pl
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub check_address
{
  $logger->trace('call MemeWebUtils->check_address') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($address, $email_contact) = @_;
  my ($status);

  if (!$address) {
    $self->whine(
      "You must include a return e-mail address to receive your results.<br />",
      "Please go back and include one."
    );
  } else {
    $status = valid_address($address);
    if ($status == 0) {
      $self->whine(
        "There is an error in your return email address:<br />",
        "&nbsp;&nbsp;&nbsp;&nbsp; <tt>$address</tt><br />",
        "It is possible that your email address is correct, in which case",
        "the problem may be that your host is behind a firewall and",
        "is consequently not found by the nslookup routines.  Consult with",
        "your systems people to see if you have an nslookup-visible email",
        "address.  If none is available, please send email to <br />",
        "&nbsp;&nbsp;&nbsp;&nbsp; <tt>$email_contact</tt> <br />",
        "mentioning this problem."
      );
    }
  }
} # check_address

#
# Check to see whether the p-value threshold is valid.
#
# Used in
# fimo.pl
# 
sub check_pvalue
{
  $logger->trace('call MemeWebUtils->check_pvalue') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($pvalue) = @_;

  if (! is_numeric($pvalue) || $pvalue < 0 || $pvalue > 1.0) {
    $self->whine(
      "You must use a number between 0 and 1 in the output threshold p-value field."
    );
  }
} # check_pvalue

#
# Check whether a DNA-only option is valid
#
# Used in
# fimo.pl
#
sub check_dna_only
{
  $logger->trace('call MemeWebUtils->check_dna_only') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($alphabet, $option, $text) = @_;

  if ($alphabet ne "DNA") {
    $self->whine("You may only use the '$option' $text.");
  }
} # check_dna_only

#
# Check to see whether the description is valid.
# 
# Used in
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub check_description
{
  $logger->trace('call MemeWebUtils->check_description') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($description) = @_;
  $description = "" unless defined $description;
  if ($description =~ /[^\w _\-\(\)]/) {
    $self->whine(
      "You may only use letters, numbers, underlines, dashes, blanks and parentheses",
      "in the description field."
    );
  }
} # check_description

#
# I wanted to use the unlink1 used by File::Temp
# but for some reason they don't allow it to be
# exported. It shouldn't matter that I don't
# get the stat compare before the delete though
# as there shouldn't be any doubt that they
# refer to the same file.
#
sub unlink1 {
  $logger->trace('call MemeWebUtils->unlink1') if $logger;
  my ($fh, $filename) = @_;
  close($fh);
  unlink($filename);
}

#
# get sequence data
#
# Object method.
#
# Get data from a textbox or uploaded file and convert
# it to FASTA format.
#
# Accepts optional parameters:
# MAX         integer     the maximum count of characters in the sequences
# SHUFFLE     boolean     should the sequences be shuffled
# PURGE       number      when set run purge with this score
# DUST        number      when set run dust with this cutoff
# NAME        string      the name of the dataset to be reported in messages
# ALPHA       DNA|PROTEIN test the alphabet of the dataset
# PASTE       boolean     when set select the pasted if true
#
# Returns fasta_data, alphabet, nseqs, min, max, ave, total
# Uses object global $self->{PROGRAM} and $self->{BIN_DIR}.
#
# Used by
# Utils
# meme.pl
# glam2.pl
# dreme.pl
#
sub get_sequence_data {
  $logger->trace('call MemeWebUtils->get_sequence_data') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my $BIN_DIR = $self->{BIN_DIR};
  my (
    $data,            # data from textbox, if defined
    $file,            # open filehandle from CGI.pm, if defined
    %options
  ) = @_;

  my $dataset_name = (defined($options{NAME}) ? $options{NAME} : 'input dataset');
  my $use_pasted = $options{PASTE};
  my $is_actual_seqs = 0;
  # return values
  my $fasta_data = ""; # sequence data in FASTA format
  my $alphabet = "UNRECOGNIZED"; # sequence alphabet DNA/PROTEIN
  my $nseqs = 0; # number of sequences
  my $min = 0;   # minimum length
  my $max = 0;   # maximum length
  my $ave = 0;   # average length
  my $total = 0; # total length
  
  # temporary files
  my ($raw_tmp, $raw_nam);
  my ($cooked_tmp, $cooked_nam);
  my ($shuffled_tmp, $shuffled_nam);
  my ($purged_tmp, $purged_nam);
  my ($dusted_tmp, $dusted_nam);
  # fasta seqs (reference to one of the above)
  my ($fasta_tmp, $fasta_nam);

  # other vars
  my ($status, $errors, $size, $has_problems);

  if (defined($use_pasted)) {
    if ($use_pasted) {
      if (!$data) {
        $self->whine(
          "No actual sequences were entered for the $dataset_name but $PROGRAM ",
          "requires the $dataset_name to run.<br />",
          "If you still wish to submit a query, please go back and either select ",
          "a sequence file to upload or paste the actual sequences for the ",
          "$dataset_name."
        );
        goto DONE;
      }
    } else {
      if (!$file) {
        $self->whine(
          "No sequence file was entered for the $dataset_name but $PROGRAM ",
          "requires the $dataset_name to run.<br />",
          "If you still wish to submit a query, please go back and either select ",
          "a sequence file to upload or paste the actual sequences for the ",
          "$dataset_name."
        );
        goto DONE;
      }
      if ($file && (-s $file) == 0) {
        $self->whine(
          "The sequence file entered for the $dataset_name is empty ",
          "but $PROGRAM requires at least 1 sequence to run.<br />",
          "If you still wish to submit a query, please go back and select ",
          "a non-empty sequence file to upload or paste actual sequences for ",
          "the $dataset_name."
        );
        goto DONE;
      }
    }
    $is_actual_seqs = $use_pasted;
  } else {
    # check that sequence data was provided
    if (!$file && !$data) {
      $self->whine(
        "No data was entered for the $dataset_name but $PROGRAM requires the ",
        "$dataset_name to run.<br />",
        "If you still wish to submit a query, please go back and either select ",
        "a sequence file to upload or paste the actual sequences for the ",
        "$dataset_name."
      );
      goto DONE;
    }

    # don't allow both datafile and data
    if ($file && $data) {
      $self->whine(
        "Both the sequence file and actual sequences were entered for the ",
        "$dataset_name.<br />",
        "If you still wish to submit a query, please go back and either clear",
        "the file selection or erase the content of the actual sequences field."
      );
      goto DONE;
    }

    # don't allow empty sequence files
    if ($file && (-s $file) == 0) {
      $self->whine(
        "The sequence file entered for the $dataset_name is empty ",
        "but $PROGRAM requires at least 1 sequence to run.<br />",
        "If you still wish to submit a query, please go back and select ",
        "a non-empty sequence file to upload or paste actual sequences for ",
        "the $dataset_name."
      );
      goto DONE;
    }
    $is_actual_seqs = !$file;
  }

  # print raw sequences to a file
  ($raw_tmp, $raw_nam) = tempfile('raw_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);
  $logger->debug('get_sequence_data - $raw_tmp: ' . $raw_nam) if $logger;
  if ($is_actual_seqs) {
    # convert to UNIX EOL
    $data =~ s/\r\n/\n/g;        # Windows -> UNIX eol
    $data =~ s/\r/\n/g;            # MacOS -> UNIX eol
    # print to file
    print $raw_tmp $data;
  } else {
    # copy the file changing to UNIX EOL
    my $line;
    while ($line = <$file>) {
     $line =~ s/\r\n/\n/g;
     $line =~ s/\r/\n/g;
     print $raw_tmp $line, "\n";
   }
  }

  # run GETSIZE on the raw sequences and see if it reports any errors
  $status = &invoke(BIN => $BIN_DIR, PROG => 'getsize', ARGS => [$raw_nam], 
    OUT_VAR => \$size, ERR_VAR => \$errors, LOGGER => $logger);
  $logger->debug('get_sequence_data - GETSIZE on raw seqs: ' . $size) if $logger;

  if ($errors || $status) {
    # maybe the file is not FASTA format? Attempt to convert it using READSEQ
    ($cooked_tmp, $cooked_nam) = tempfile('cooked_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);
    $status = &invoke(BIN => $BIN_DIR, PROG => 'readseq', ARGS => 
      ['-a', '-f=8', $raw_nam], OUT_FILE => $cooked_nam, ERR_VAR => \$errors);

    # check for errors
    if ($status) {
      $self->whine(
        "The sequences submitted for the $dataset_name could not be read as ",
        "FASTA format and automatic conversion using the READSEQ program ",
        "failed.<br />",
        "READSEQ returned: <pre>$errors</pre><br />",
        "If you still wish to submit a query, please go back and select ",
        "a FASTA formatted sequence file to upload or paste actual ",
        "sequences in FASTA format for the $dataset_name."
      );
      goto DONE;
    }

    #run GETSIZE on the converted sequences
    #(note -nd means do not print warnings about duplicate sequences)
    $status = &invoke(BIN => $BIN_DIR, PROG => 'getsize', ARGS => 
        ['-nd', $cooked_nam], OUT_VAR => \$size, 
        ERR_VAR => \$errors);
    $logger->debug('get_sequence_data - GETSIZE on cooked seqs: ' . $size) if $logger;
  
    # check for errors
    if ($errors || $status) {
      $self->whine(
        "The sequences submitted for the $dataset_name were converted to ",
        "FASTA format but in the process of detecting the alphabet and size ",
        "the following errors in your dataset were detected:<br /><pre>$errors</pre>",
        "<br />Make sure all your sequences are in the same format since READSEQ",
        "assumes that all sequences are in the same format as the first sequence."
      );
      goto DONE;
    }
    # use converted dataset
    ($fasta_tmp, $fasta_nam) = ($cooked_tmp, $cooked_nam);
    &unlink1($raw_tmp, $raw_nam); # this should be done when perl exits as well
  } else {
    # use original dataset
    ($fasta_tmp, $fasta_nam) = ($raw_tmp, $raw_nam);
  }
  # just so we don't close a file twice later
  $raw_tmp = undef; $raw_nam = undef; $cooked_tmp = undef; $cooked_nam = undef;

  # extract out the sequence stats from getsize's response
  my $letters;
  ($nseqs, $min, $max, $ave, $total, $letters) = split (' ', $size);

  $has_problems = 0;

  # check for problem reading the dataset
  if ($nseqs <= 0) {
    $self->whine(
      "The sequences submitted for the $dataset_name appear to be ",
      "in a format that $PROGRAM does not recognize.<br /> ",
      "Please check to be sure that your data is ",
      "<a href=\"../help_format.html\">formatted</a> properly."
    );
    $has_problems = 1;
  }
  # check for bad sequences
  if ($nseqs > 0 && $min == 0) {
    $self->whine(
      "The sequences submitted for the $dataset_name appear to ",
      "contain one or more zero-length sequences.<br /> ",
      "Please check to be sure that your data is ",
      "<a href=\"../help_format.html\"> formatted</a> properly."
    );
    $has_problems = 1;
  }
  # Make sure there isn't too much data.
  if (defined($options{MAX}) && $total > $options{MAX}) {
    $self->whine(
      "The sequences submitted for the $dataset_name contain ",
      &commify($total) . " characters but $PROGRAM can only accept ",
      &commify($options{MAX}) . " characters.<br/>",
      "Please submit a smaller dataset."
    );
    $has_problems = 1;
  }

  # prepare sequences and alphabet and do filtering if requested

  # calculate the alphabet
  my $bad_symbols;
  ($alphabet, $bad_symbols) = get_alph($letters); # here alphabet is "DNA" or "PROTEIN" unless error in input

  # check the alphabet
  if ($bad_symbols) {
    # load in the fasta sequences
    $data = do {local $/; <$fasta_tmp>};

    my @bad_lines = $self->find_bad_sequence_data($bad_symbols, $data);
    $self->whine(
      "The sequences submitted for the $dataset_name contained the ",
      "following unrecognized letters: $bad_symbols.<br/>",
      "The unrecognized letters occurred in the following locations:",
      @bad_lines,
      "<br/>",
      "Please convert your sequences to one of the recognized sequence",
      "<a href=\"../help_alphabet.html\">alphabets</a>."
    );
    $has_problems = 1;
  }

  if (defined($options{ALPHA}) && $alphabet ne $options{ALPHA}) {
    $self->whine(
      "The sequences submitted for the $dataset_name used the $alphabet ",
      "alphabet but $PROGRAM expected the " . $options{ALPHA} . " alphabet."
    );
    $has_problems = 1;
  }

  goto DONE if $has_problems;

  # shuffle sequences if requested
  if ($options{SHUFFLE}) {
    ($shuffled_tmp, $shuffled_nam) = tempfile('shuffled_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);

    $status = &invoke(BIN => $BIN_DIR, PROG => 'fasta-shuffle-letters', ARGS => ['-tod'], 
        IN_FILE => $fasta_nam, OUT_FILE => $shuffled_nam, ERR_VAR => \$errors);

    if ($errors || $status) {
      $self->whine(
        "When shuffling the $dataset_name, the following errors ",
        "resulted:<br /><pre>$errors</pre>",
        "<br />Please check your sequences."
      );
      goto DONE;
    }
    &unlink1($fasta_tmp, $fasta_nam); # this should be done when perl exits as well
    ($fasta_tmp, $fasta_nam) = ($shuffled_tmp, $shuffled_nam);
    $shuffled_tmp = undef;
    $shuffled_nam = undef;
  }

  # purge: remove overly similar sequences
  if (defined($options{PURGE})) {
    my @purge_args = ($fasta_nam, $options{PURGE});
    push(@purge_args, '-n') if ($alphabet eq "DNA");

    # purge is annoying... I can not direct it to write to a temporary file, 
    # instead it always writes to a file called <input file>.<purge score>
    my $purge_out = $fasta_nam . '.' . $options{PURGE};
    log_and_die("Purge output file already exists!") if (-e $purge_out);
    $status = &invoke(BIN => $BIN_DIR, PROG => 'purge', ARGS => \@purge_args, ERR_VAR => \$errors); 
    # check for errors
    if ($errors || $status) {
      $self->whine("When calling purge on the $dataset_name, the following ",
        "errors resulted:<br /><pre>$errors</pre>");
      unlink($purge_out) if (-e $purge_out);
      goto DONE;
    }
    # make things neater by copying the file created by purge 
    # to a temporary file and deleteing the original
    ($purged_tmp, $purged_nam) = tempfile('purged_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);
    &copy($purge_out, $purged_tmp);
    unlink($purge_out);

    &unlink1($fasta_tmp, $fasta_nam); # this should be done when perl exits as well
    ($fasta_tmp, $fasta_nam) = ($purged_tmp, $purged_nam);
    $purged_tmp = undef;
    $purged_nam = undef;
  }

  # dust: removes low information content regions
  if (defined($options{DUST})) {
    unless ($alphabet eq "DNA") {
      $self->whine ("Could not call dust on the $dataset_name, ",
        "dust is only good for DNA");
      goto DONE;
    }
    # make a temp output file
    ($dusted_tmp, $dusted_nam) = tempfile('dusted_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);
    # run dust
    $status = &invoke(BIN => $BIN_DIR, PROG => 'dust', ARGS => 
      [$fasta_nam, $options{DUST}], OUT_FILE => $dusted_nam, ERR_VAR => \$errors);
    # check for errors
    if ($errors || $status) {
      $self->whine("When calling dust on the $dataset_name, the following ",
        "errors resulted:<br /><pre>$errors</pre>");
      goto DONE;
    }

    &unlink1($fasta_tmp, $fasta_nam); # this should be done when perl exits as well
    ($fasta_tmp, $fasta_nam) = ($dusted_tmp, $dusted_nam);
    $dusted_tmp = undef;
    $dusted_nam = undef;
  }

  seek($fasta_tmp, SEEK_SET, 0); # I'm not sure if this is needed
  $fasta_data = do {local $/; <$fasta_tmp>};

DONE:
  # clean up
  &unlink1($raw_tmp, $raw_nam) if $raw_tmp;
  &unlink1($cooked_tmp, $cooked_nam) if $cooked_tmp;
  &unlink1($shuffled_tmp, $shuffled_nam) if $shuffled_tmp;
  &unlink1($purged_tmp, $purged_nam) if $purged_tmp;
  &unlink1($dusted_tmp, $dusted_nam) if $dusted_tmp;
  &unlink1($fasta_tmp, $fasta_nam) if $fasta_tmp;
  # return sequences and attributes
  return($fasta_data, $alphabet, $nseqs, $min, $max, $ave, $total);
} # get sequence data

#
# get fasta data
#
# Object method.
#
# Get data from a textbox or uploaded file in FASTA format
#
# Returns: data, undef, nseqs, min, max, ave, total
# Uses object globals $self->{PROGRAM}, $self->{BIN_DIR}.
#
# meme-seq.pl
#
sub get_fasta_data {
  $logger->trace('call MemeWebUtils->get_fasta_data') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my $BIN_DIR = $self->{BIN_DIR};
  my (
    $data,   # data from textbox, if defined
    $file,   # open filehandle from CGI.pm, if defined
    $expected_alphabet, # if defined, string of only accepted sequence letters
    $expected_alphabet_name # if defined name of alphabet
  ) = @_;
  # set default
  $expected_alphabet_name = "expected" unless $expected_alphabet_name;
  # return values
  my $nseqs = 0; # number of sequences
  my $min = 0;   # minimum length
  my $max = 0;   # maximum length
  my $ave = 0;   # average length
  my $total = 0; # total length

  # check that sequence data was provided
  if (!$file && !$data) {
    $self->whine(
      "You haven't entered any sequence data. <br>",
      "If you still wish to submit a query, please go back and enter the",
      "name of a sequence file or the actual sequences."
    );
    goto DONE;
  }

  # don't allow both datafile and data
  if ($file && $data) {
    $self->whine(
      "You may not enter <i>both</i> the name of a sequence file and sequences.<br>",
      "If you still wish to submit a query, please go back and erase either",
      "what you have written in the <i>name of a file</i> field or",
      "in the <i>actual sequences</i> field."
    );
    goto DONE;
  }

  if ($file && (-s $file) == 0) {
    $self->whine("Your sequence file is empty.");
    goto DONE;
  }

  #
  # get information on sequences
  #
  my($status, $getsize_out, $getsize_err, $letters);

  # slurp the uploaded file into a scalar if there was no textbox data
  $data = do {local $/; <$file>} unless ($data);

  # convert to UNIX EOL
  $data =~ s/\r\n/\n/g;        # Windows -> UNIX eol
  $data =~ s/\r/\n/g;            # MacOS -> UNIX eol

  # write the data to a temp file
  my ($seq_tmp, $seq_nam) = tempfile('raw_seqs_XXXXXXXXXX', DIR => $tmpdir, UNLINK => 1);
  print $seq_tmp $data;

  # Run the 'getsize' program to get information on the raw sequence data.
  $status = &invoke(BIN => $BIN_DIR, PROG => 'getsize', ARGS => [$seq_nam], 
    OUT_VAR => \$getsize_out, ERR_VAR => \$getsize_err);

  # close and remove the temp file
  &unlink1($seq_tmp, $seq_nam);

  # handle errors
  if ($getsize_err || $status) {
    $self->whine(
      "Check your FASTA sequences. ",
      "The following errors in your dataset were detected:<br>",
      '<div style="font-family:monspace;">',
      map($_ . '<br>', split(/\n/, $getsize_err)),
      "</div>"
    );
    goto DONE;
  }

  # extract the information returned by getsize
  ($nseqs, $min, $max, $ave, $total, $letters) = split (' ', $getsize_out);

  #
  # check the information we gathered
  #

  # if we have an expected alphabet, whine if a letter isn't in it
  # in any sequence (case independent)
  if ($expected_alphabet && $letters !~ m/^[$expected_alphabet]+$/i) {
    $self->whine(
      "Check your FASTA sequences: they contain at least one<br>".
      "character not in the $expected_alphabet_name alphabet<br>".
      "(allowed characters: `$expected_alphabet')."
    );
    goto DONE;
  }
  # check for bad sequences
  if ($nseqs > 0 && $min == 0) {
    $self->whine(
      "Your dataset appears to contain one or more zero-length sequences.",
      "<br /> Please check to be sure that your data is",
      "<a href=\"../doc/fasta-format.html\">formatted</a> properly."
    );
    goto DONE;
  }
  # Make sure the data was in FASTA format or got converted correctly.
  if ($self->{NERRORS} == 0 && $nseqs == 0 ) {
    $self->whine(
      "$PROGRAM was unable to read your sequence data.<br />",
      "Please check to be sure that your sequence data is in correct",
      '<a href=\"../doc/fasta-format.html\">FASTA</a> format.<br />',
      "If you are still having trouble, you can try to convert your data to",
      '<a href=\"../doc/fasta-format.html\">',
      "FASTA format</a> using the",
      '<a href=\"../doc/readseq.html\">',
      "ReadSeq</a> program and then resubmit it."
    );
  }

DONE:
  my $alphabet = undef; # apparently not used!
  return($data, $alphabet, $nseqs, $min, $max, $ave, $total);
} # get fasta data


#
# Get actual name of the database to search along with
# other information about it.
#
# Uses global $PROGRAM
#
# Used by
# fimo.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub get_db_name {
  $logger->trace('call MemeWebUtils->get_db_name') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my (
    $motif_alphabet,    # motif alphabet
    $upload_db_name,    # name of uploaded db file
    $max_upload_size,   # maximum size of the uploaded file
    $translate_dna,    # scan DNA with protein motifs
    $short_dna_only,    # don't allow DNA search of long sequences
    $database_index,    # index of databases
    $database_id        # id of database
  ) = @_;
  my (
    $db,         # full name of db file
    $uploaded_data,    # contents of uploaded file
    $db_alphabet,    # local/PROTEIN/DNA
    $target_alphabet,    # alphabet of db
    $db_root,         # db file name w/o extension
    $prot_db,         # protein version exists
    $nuc_db,         # DNA version exists
    $short_seqs,    # sequences are all "short"
    $db_menu_name,    # db name for menus
    $db_long_name,    # db name for documentation
    $db_descr,         # text for verification message
    $prom_start,    # start of promoter (gomo)
    $prom_end,        # end of promoter (gomo)
    @extra_dbs        # extra databases (gomo)
  );    # return values

  # get available versions of database
  if ($upload_db_name eq "") { # local database
    $db_alphabet = "local";
    #attempt to load the entry
    eval {
      ($db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $prom_start, $prom_end, @extra_dbs) =
          load_entry($database_index, $database_id);
    };
    if ($@) {
      #loading most likely failed because they didn't select a database
      $self->whine(
        "You must specify a supported database or a database to upload.<br />",
        "Please go back and specify one."
      );
    }
    $db_descr = "";
  } else { # uploaded database
    $db_root = $upload_db_name;
    # check that name doesn't contain quotes
    if ($db_root =~ /['`"]/) {
      $self->whine("Uploaded database file name may not contain quote characters.");
    }
    my ($nseqs, $min, $max, $ave, $total);
    ($uploaded_data, $db_alphabet, $nseqs, $min, $max, $ave, $total) =
        $self->get_sequence_data("", $upload_db_name, $max_upload_size, 0);
    $prot_db = ($db_alphabet eq "PROTEIN") ? 1 : 0;
    $nuc_db = ($db_alphabet eq "DNA") ? 1 : 0;
    $short_seqs = 1;                # allow search 
    $db_descr = "  <li>Database type: <b>$db_alphabet</b></li>\n".
    "<li>Database sequences: <b>$nseqs</b></li>\n".
    "<li>Database size: <b>$total</b></li>";
    #does this mean there are no extra dbs? because the user uploaded this db? I'm not certain
    @extra_dbs = ();
  }
  my $i;
  my @extra_files = ();
  while (@extra_dbs) {
    push(@extra_files, shift(@extra_dbs) . ".na");
    shift(@extra_dbs); #skip description
  }
  # get desired sequence alphabet
  $target_alphabet = ($motif_alphabet eq "PROTEIN" && !$translate_dna) ? "PROTEIN" : "DNA";

  # set actual name of database to search
  if ($db_alphabet eq "local") {
    $db = $db_root;
    if ($target_alphabet eq "DNA" && $nuc_db) {
      if ($short_dna_only && !($short_seqs)) {
        $self->whine("The DNA sequences in database '$db_menu_name' are too long for searching by $PROGRAM.");
      } else {
        $db .= ".na";
      }
    } elsif ($target_alphabet eq "PROTEIN" && $prot_db) {
      $db .= ".aa";
    } elsif (defined $db_root) {
      $self->whine("There is no $target_alphabet version of database '$db_menu_name'.");
    }
  } else {
    if ($target_alphabet eq $db_alphabet) {
      $db = "uploaded_db";
    } else {
      $self->whine("The uploaded database '$db_root' seems to be $db_alphabet but should be $target_alphabet.");
    }
  }

  return($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db,
    $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr, \@extra_files);
} # get_db_name

#
# Upload the file containing motifs and return data in a scalar.
# Converts to UNIX EOL.  Checks that the file is not empty.
# 
# Returns $motif_data.
#
# Used by
# fimo.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub upload_motif_file {
  $logger->trace('call MemeWebUtils->upload_motif_file') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my (
    $motif_dualval_file,    # (raw) uploaded motif file and filename in a dualval
    $inline_motifs,         # contains inline motifs if defined (input)
    $choose_inline          # force choice of one if other is avaliable
  ) = @_;
  my $motif_file = openhandle($motif_dualval_file);
  my $motif_data = $inline_motifs;

  if (defined $choose_inline) {
      if ($choose_inline) {
        unless ($inline_motifs) {
          $self->whine("The inline motifs/alignment were not provided so you ".
            "can not choose them. This shouldn't happen when the form is ".
            "designed properly so unless you're using a non-standard form ".
            "this is probably a bug.");
        }
      } else {
        unless ($motif_file) {
          $self->whine("No motif/alignment file was uploaded.");
        }
        $motif_data = "";
      }
  }
  # slurp the uploaded file into a scalar if there was no inline data
  $motif_data = do {local $/; <$motif_file>} unless ($motif_data);

  # convert to UNIX EOL
  $motif_data =~ s/\r\n/\n/g;        # Windows -> UNIX eol
  $motif_data =~ s/\r/\n/g;        # MacOS -> UNIX eol

  #
  # check that motif file was found and not empty
  #
  unless (defined $motif_data && $motif_data =~ /\S+/) {
    $self->whine("Your motif/alignment file ".
      "was empty.  Please check it and retry.");
    $motif_data = "";
  }

  return($motif_data);
} # upload_motif_file

#
# Check the validity of a MEME file.
#
# Returns $motifs_found, $motif_alphabet, $nmatrices and $total_cols
#
# Accepts options:
# ALPHA   DNA|PROTEIN   Alphabet required by motifs.
# NAME    string        Name of the dataset for error messages.
# 
# Used by
# fimo.pl
# gomo.pl
# mast.pl
# mcast.pl
# tomtom.pl
#
sub check_meme_motifs {
  $logger->trace('call MemeWebUtils->check_meme_motifs') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my $BIN_DIR = $self->{BIN_DIR};
  my (
    $motif_data,          # string containing motifs
    %options
  ) = @_;
  my ($motifs_found, $motif_alphabet, $total_cols, $nmatrices, $version);
  my ($alphabet);
  my $version_hdr_re = 'MEME\s+version\s+(\S+)';
  my $alph_hdr_re = 'ALPHABET=\s*(\S+)';
  my $strand_hdr_re = 'strands:';
  my $motif_hdr_re = 'letter-probability\s+matrix:\s.*\sw\s*=\s+(\d+)';
  my $dataset_name = (defined($options{NAME}) ? $options{NAME} : 'input motifs');

  # process meme-io type file
  # first try parsing as XML
  ($motifs_found, $motif_alphabet, $nmatrices, $total_cols) =
      $self->check_meme_motifs_as_xml($motif_data);

  if (!$motifs_found) {
    # second try parsing as HTML
    ($motifs_found, $motif_alphabet, $nmatrices, $total_cols) =
        $self->check_meme_motifs_as_html($motif_data);
  }

  # otherwise, try parsing as HTML or text
  if (!$motifs_found) {
    my @widths;
    # check that MEME version given
    ($version) = ($motif_data =~ /$version_hdr_re/);
    $self->whine("Can't find MEME version in $dataset_name file.") unless ($version);
    # get the alphabet
    ($alphabet) = ($motif_data =~ /$alph_hdr_re/);
    chop($_ = `$BIN_DIR/alphtype $alphabet 2>&1`);
    $motif_alphabet = $_;                     # DNA/PROTEIN
    # check if strands given for DNA
    if ($motif_alphabet eq 'DNA') {
      my ($strands) = ($motif_data =~ /$strand_hdr_re/);
      $self->whine("Can't find 'strands:' in $dataset_name file.") unless ($strands);
    }
    # get the number of motifs, total columns
    $nmatrices = @widths = ($motif_data =~ /$motif_hdr_re/g);
    $total_cols = eval(join "+", @widths);
    $motifs_found = ($nmatrices > 0);
  }
  if ($motif_alphabet ne "PROTEIN" and $motif_alphabet ne "DNA") {
    $self->whine(
      "The motifs submitted for the $dataset_name used an unknown alphabet."
    );
    $motifs_found = 0;
    $motif_alphabet = 'UNKNOWN';
  } elsif (defined($options{ALPHA}) && $motif_alphabet ne $options{ALPHA}) {
    $self->whine(
      "The motifs submitted for the $dataset_name used the $motif_alphabet ",
      "alphabet but $PROGRAM expected the " . $options{ALPHA} . " alphabet."
    );
    $motifs_found = 0;
  }

  return($motifs_found, $motif_alphabet, $nmatrices, $total_cols);

} # check_meme_motifs 

#
# Private Method
#
# Check the validity of a MEME XML file.
#
# Returns $motifs_found, $motif_alphabet, $num_motifs and $total_motif_columns
# 
# Used by
# Utils
#
sub check_meme_motifs_as_xml {
  $logger->trace('call MemeWebUtils->check_meme_motifs_as_xml') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $BIN_DIR = $self->{BIN_DIR};

  my $motif_data = pop @_;
  my ($motifs_found, $motif_alphabet, $num_motifs, $total_motif_columns);
  my $ref = eval { XMLin($motif_data, KeepRoot => 1) }; # Try parsing as XML

  if ($@) {
    # Parsing as XML failed.
    $motifs_found = 0;
  }
  else {
    # Parsing as XML succeded.
    if ($ref->{MEME}->{training_set}->{alphabet}) {
      # Probably a MEME XML file
      my $alphabet = $ref->{'MEME'}->{'training_set'}->{'alphabet'}; 
      $motifs_found = 1;
      my $letters = $alphabet->{'letter'};
      my @letters;
      while (my ($key, $letter) = each %$letters) {
        push @letters, "$letter->{'symbol'}";
      }
      @letters = sort @letters;
      my $alphabet_string = join '', @letters;
      chop($motif_alphabet = `$BIN_DIR/alphtype $alphabet_string 2>&1`);
      
      if ($motif_alphabet ne "PROTEIN" and $motif_alphabet ne "DNA") {
        $self->whine("Untranslatable alphabet string in XML: $alphabet_string");
        $motif_alphabet = 'UNKNOWN';
      }

      my $model = $ref->{'MEME'}->{'model'};
      $num_motifs = $model->{'nmotifs'};
      my $motifs = $ref->{'MEME'}->{'motifs'}->{'motif'};
      $total_motif_columns = 0;
      while (my ($key, $motif) = each %$motifs) {
        $total_motif_columns += $motif->{'width'};
      }
    } elsif ($ref->{dreme}->{motifs}) {
      # Probably a DREME XML file
      $motifs_found = 1;
      my $type = $ref->{dreme}->{model}->{background}->{type};
      if ($type ne 'dna') {
        $self->whine("Unhandled alphabet type in DREME XML: $type");
      }
      $motif_alphabet = 'DNA';
      my $motifs = $ref->{dreme}->{motifs}->{motif};
      $total_motif_columns = 0;
      $num_motifs = 0;
      while (my ($key, $motif) = each %$motifs) {
        $total_motif_columns += $motif->{length};
        $num_motifs++;
      }
    } else {
      # Unknown XML file
      $motifs_found = 0;
    }
  }
  return($motifs_found, $motif_alphabet, $num_motifs, $total_motif_columns);
} # check_meme_motifs_as_xml

#
# Private Method
#
# Check the validity of a MEME HTML file.
#
# Returns $motifs_found, $motif_alphabet, $num_motifs and $total_motif_columns
# 
# Used by
# Utils
#
sub check_meme_motifs_as_html {
  $logger->trace('call MemeWebUtils->check_meme_motifs_as_html') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $BIN_DIR = $self->{BIN_DIR};
  my $motif_data = pop @_;
  my ($motifs_found, $motif_alphabet, $num_motifs, $total_motif_columns);
  my $motif_hdr_re = 'letter-probability matrix:\s+alength=\s+\d+\s+w=\s+(\d+)';

  my $p = HTML::PullParser->new(doc => \$motif_data, start => 'attr', report_tags => qw(input));
  my $token = $p->get_token;
  while ($token) {
    my $hashref = @$token[0]; 
    my %attrs = %$hashref;
    if (exists($attrs{"type"}) && $attrs{"type"} =~ m/^hidden$/i && exists($attrs{"name"}) && exists($attrs{"value"})) {
      my $name = $attrs{"name"};
      my $value = $attrs{"value"};
      if ($name =~ m/^nmotifs$/i) {
        $num_motifs = int($value);
        if ($num_motifs <= 0) {
          $self->whine("The number of motifs was smaller than expected: $value");
        }
      } elsif ($name =~ m/^alphabet$/i) {
        chop($motif_alphabet = `$BIN_DIR/alphtype $value 2>&1`);
        if ($motif_alphabet ne "PROTEIN" and $motif_alphabet ne "DNA") {
          $self->whine("Untranslatable alphabet string in HTML: $value");
          $motif_alphabet = 'UNKNOWN';
        }
      } elsif ($name =~ m/^pspm\d+$/i) {
        my $width = int($value =~ m/letter-probability matrix:\s+alength=\s+\d+\s+w=\s+(\d+)/i);
        if (defined $total_motif_columns) {
          if ($width != $total_motif_columns) {
            $self->whine("$name has a different width ($width) to the expected ($total_motif_columns)");
          }
        } else {
          $total_motif_columns = $width;
        }
      } 
    }
    $token = $p->get_token;
  }
  if (defined $num_motifs && defined $motif_alphabet && defined $total_motif_columns) {
    $motifs_found = 1;
  } else {
    $motifs_found = 0;
  }

  return($motifs_found, $motif_alphabet, $num_motifs, $total_motif_columns);
} # check_meme_motifs_as_html

#
# Private Method
#
# send a verification message
#
# Uses global $PROGRAM
#
# Used by
# Utils
#
sub send_verification
{
  $logger->trace('call MemeWebUtils->send_verification') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my (
    $address,            # to
    $from,            # to
    $jobid,            # from opal
    $message             # message
  ) = @_;

  my $sendmail = "/usr/sbin/sendmail -t";

  my $content = "To: " . $address . "\n";
  $content .= "From: $from\n";
  $content .= "Reply-to: $from\n";
  $content .= "Subject: $PROGRAM Submission Information (job $jobid)\n";
  $content .= "Content-type: text/html\n\n";
  $content .= "<hr>This is an auto-generated response to your job submission.<br />\n\n";
  $content .= $message;

  open(SENDMAIL, "|$sendmail") or $self->whine("Can't open sendmail.<br />");
  print SENDMAIL $content;
  close SENDMAIL;
} # send_verification


#
# Get a required parameter and die if it's not submitted
#
sub param_require {
  $logger->trace('call MemeWebUtils->param_require') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  my $data = $q->param($name);
  unless (defined $data) {
    $logger->error("Required field '$name' was not passed to the CGI script.") if $logger;
    die("Required field '$name' was not passed to the CGI script.\n");
  } else {
    $logger->debug("Required field '$name' exists.") if $logger;
  }
  return $data;
}

#
# Get a choice parameter
# The parameter must take one of the choices or the script dies.
#
sub param_choice {
  $logger->trace('call MemeWebUtils->param_choice') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name, @choices) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  my $data = $self->param_require($q, $name);
  foreach (@choices) {
    if ($data eq $_) {
      $logger->debug("Choice field '$name' = '$data'") if $logger;
      return $_;
    }
  }
  log_and_die("Choice parameter did not return one of the valid choices");
}

#
# Get a boolean parameter, undefined counts as false
#
sub param_bool {
  $logger->trace('call MemeWebUtils->param_bool') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  my $data = $q->param($name);
  if (not defined $data) {
    $logger->debug("Bool field '$name' = false") if $logger;
    return 0;
  }
  if ($data eq '1') {
    $logger->debug("Bool field '$name' = true") if $logger;
    return 1;
  } else {
    $logger->debug("Bool field '$name' = false") if $logger;
    return 0;
  }
}

#
# Get a num parameter, empty or whitespace returns the default
#
sub param_num {
  $logger->trace('call MemeWebUtils->param_num') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name, $desc, $min, $max, $default) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  my $data = $self->param_require($q, $name);
  if ($data !~ m/^\s*$/) {
    my $num = getnum($data);
    if (defined($num)) {
      if (defined($min) && $num <= $min) {
        $self->whine("The $desc ($num) must be &gt; $min.");
      } elsif (defined($max) && $num > $max) {
        $self->whine("The $desc ($num) must be &le; $max.");
      } else {
        $logger->debug("Num field '$name' = $num") if $logger;
        return $num;
      }
    } else {
      $self->whine("The $desc must be a number.");
    }
  }
  $logger->debug("Num field '$name' = $default (default)") if $logger;
  return $default;
}

#
# Get a int parameter, empty or whitespace returns the default
#
sub param_int {
  $logger->trace('call MemeWebUtils->param_int') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name, $desc, $min, $max, $default) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  my $data = $self->param_require($q, $name);
  if ($data !~ m/^\s*$/) {
    if ($data =~ m/^\s*([+-]?\d+)\s*$/) {
      my $num = $1;
      if (defined($min) && $num < $min) {
        $self->whine("The $desc ($num) must be &ge; $min.");
      } elsif (defined($max) && $num > $max) {
        $self->whine("The $desc ($num) must be &le; $max.");
      } else {
        $logger->debug("Int field '$name' = $num") if $logger;
        return $num;
      }
    } else {
      $self->whine("The $desc must be a number.");
    }
  }
  $logger->debug("Int field '$name' = $default (default)") if $logger;
  return $default;
}

#
# Get a series of parameters which make a background
#
sub param_dna_bg {
  $logger->trace('call MemeWebUtils->param_dna_bg') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $name, $desc) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  my %bg = (
    A => $self->param_num($q, $name . 'A', $desc . ' A', 0, undef, 0.25), 
    C => $self->param_num($q, $name . 'C', $desc . ' C', 0, undef, 0.25), 
    G => $self->param_num($q, $name . 'G', $desc . ' G', 0, undef, 0.25), 
    T => $self->param_num($q, $name . 'T', $desc . ' T', 0, undef, 0.25));
  #normalise the background
  my $bg_sum = $bg{A} + $bg{C} + $bg{G} + $bg{T};
  $bg{A} /= $bg_sum;
  $bg{C} /= $bg_sum;
  $bg{G} /= $bg_sum;
  $bg{T} /= $bg_sum;
  $bg{dna} = 1;
  $bg{source} = "web form";
  return \%bg;
}

#
# Get a filter parameter
#
sub param_filter {
  $logger->trace('call MemeWebUtils->param_filter') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, $param_all, $param_filter, $desc) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';

  return (0, [], []) if ($self->param_bool($q, $param_all));

  my @index_filter = ();
  my @name_filter = ();
  foreach my $id (split(/\s/, $self->param_require($q, $param_filter))) {
    if ($id =~ m/^\d+$/) {
      push(@index_filter, int($id));
    } elsif (is_safe_name($id)) {
      push(@name_filter, $id);
    } else {
      $self->whine("The value \"$id\" is not an allowed name ".
        "for the $desc filter. You may use alphanumeric characters with ".
        "colons ':', underscores '_', dots '.' and dashes '-' ".
        "(provided the dash isn't at the front).");
    }
  }
  my $filter_count = scalar(@index_filter) + scalar(@name_filter);
  return ($filter_count, \@index_filter, \@name_filter);
}

#
# Extracts motif database parameters from a passed CGI object
#
# Accepts the following options
# DB    the name of the parameter which contains the database index or 'upload' or 'inline'
# UP    the name of the parameter which contains the file when 'upload' is specified
# IN    the name of the parameter which contains inline motifs when 'inline' is specified
# FNAME the name of the parameter which contains the source file name of inline motifs
# ALT   the alternate filename to use when the file parameter has an invalid name
# NAME  the descriptive name of the motif database
# ALPHA the required motif alphabet DNA or PROTEIN
#
# returns db_name db_databases up_origname up_name up_data 
#         up_alphabet up_count up_cols
#
sub param_motif_db {
  $logger->trace('call MemeWebUtils->param_motif_db') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, %opts) = @_;
  my ($db_name, $db_databases, $up_origname, $up_name, $up_data);
  my ($up_alphabet, $up_count, $up_cols) = ('UNKNOWN', 0, 0);

  my $param_db = (defined($opts{DB}) ? $opts{DB} : 'motif_db');
  my $param_file = (defined($opts{UP}) ? $opts{UP} : 'motif_db_upload');
  my $param_inline = (defined($opts{IN}) ? $opts{IN} : 'inline_motifs');
  my $param_inline_name = (defined($opts{FNAME}) ? $opts{FNAME} : 'inline_name');
  my $alt = (defined($opts{ALT}) ? $opts{ALT} : 'motif_db.meme');
  my $name = (defined($opts{NAME}) ? $opts{NAME} : 'motif database');
  my $db_index = (defined($opts{DB_INDEX}) ? $opts{DB_INDEX} : $MOTIF_DB_INDEX);

  my $db = $self->param_require($q, $param_db);
  if ($db eq 'upload') {
    my $file = $q->upload($param_file);
    unless (defined($file)) {
      $self->whine('You have selected to upload the '.$name.' file'.
        'but have not selected a file to upload.');
    } else {
      $up_data = do { local $/;  <$file> };
      close($file);
      $up_origname = fileparse($q->param($param_file));
      $up_name = get_safe_name($up_origname, $alt, 0);
      # check motifs
      (undef, $up_alphabet, $up_count, $up_cols) =
      $self->check_meme_motifs($up_data, NAME => $name, ALPHA => $opts{ALPHA});
    }
  } elsif ($db eq 'inline') {
    $up_data = $q->param($param_inline);
    $up_origname = fileparse($q->param($param_inline_name));
    $up_name = get_safe_name($up_origname, $alt, 0);
    # check motifs
    (undef, $up_alphabet, $up_count, $up_cols) =
    $self->check_meme_motifs($up_data, NAME => $name, ALPHA => $opts{ALPHA});
  } else {
    if ($db =~ m/^(\d+)$/) {
      my $index = int($1);
      my @entry = load_entry($db_index, $index);
      $db_name = $entry[CatList::DB_NAME];
      $db_databases = $entry[CatList::DB_ROOT];
    } else {
      $self->whine('You must select a '.$name.'.');
    }
  }
  return ($db_name, $db_databases, $up_origname, $up_name, $up_data, 
    $up_alphabet, $up_count, $up_cols);
}

sub param_email {
  $logger->trace('call MemeWebUtils->param_email') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my ($q, %opts) = @_;
  log_and_die("Expected CGI object") unless ref($q) eq 'CGI';
  # extract options and use defaults where required
  my $param_email = (defined($opts{EMAIL}) ? $opts{EMAIL} : 'email');
  my $param_confirm = (defined($opts{CONFIRM}) ? $opts{CONFIRM} : 'email_confirm');
  my $help = (defined($opts{HELP}) ? $opts{HELP} : '@contact@');

  # get the email address
  my $email = $self->param_require($q, $param_email);
  my $email_confirm = $self->param_require($q, $param_confirm);
  if ($email ne $email_confirm) {
    $self->whine("Your email address does not match the confirmation email address.");
  } else {
    $self->check_address($email, $help);
  }
  return $email;
}

#
# make form header
# 
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_form_header {
  $logger->trace('call MemeWebUtils->make_form_header') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($program, $form_type, $optional_css, $indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;

  my $content = 
      $indent."<!DOCTYPE html>\n".
      $indent."<html>\n".
      $indent.$tab."<head>\n".
      $indent.$tab.$tab."<title>$program - $form_type </title>\n".
      $indent.$tab.$tab."<!-- provide relative path to html directory (parent of cgi-bin) -->\n".
      $indent.$tab.$tab."<script language=\"javascript\" type=\"text/javascript\">var html_path = \"../\"</script>\n".
      $indent.$tab.$tab."<script src=\"../check-submission-form.js\" type=\"text/javascript\"></script>\n".
      $indent.$tab.$tab."<script src=\"../tomtom.js\" type=\"text/javascript\"></script>\n".
      $indent.$tab.$tab."<link href=\"../doc/meme-suite.css\" rel=\"stylesheet\" type=\"text/css\" />\n".
      $indent.$tab.$tab."<link rel=\"icon\" type=\"image/png\" href=\"../images/memesuite_icon.png\">" .
      $indent.$tab.$tab."<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"../images/memesuite_icon.ico\">" .
      $indent.$tab.$tab."<!-- provide function for clearing upload field -->\n".
      $indent.$tab.$tab."<script> function clearFileInputField(tagId) { document.getElementById(tagId).innerHTML = document.getElementById(tagId).innerHTML; } </script>\n";
  if ($optional_css) {
      $content .=
          $indent.$tab.$tab."<style>\n".
          $optional_css.
          $indent.$tab.$tab."</style>\n";
  }
  $content .=
      $indent.$tab."</head>\n";

 return($content);
} # make_form_header

#
# make submission form tailer
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_submission_form_tailer {
  $logger->trace('call MemeWebUtils->make_submission_form_tailer') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content = 
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.
      "<script src=\"../template-footer.js\" type=\"text/javascript\" ></script>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.
      "<script type=\"text/javascript\" >make_footer('../');</script>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."</div>\n".
      $indent.$tab.$tab.$tab.$tab."</td>\n".
      $indent.$tab.$tab.$tab."</tr>\n".
      $indent.$tab.$tab."</table>\n".
      $indent.$tab."</body>\n".
      $indent."</html>\n";

  return($content);
} # make_submission_form_tailer

#
# make submissiom form top: program name, logo, description
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_submission_form_top {
  $logger->trace('call MemeWebUtils->make_submission_form_top') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($action, $logo, $alt, $description, $indent_lvl, $tab) = @_;
  my ($content);

  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;

  my $indent = $tab x $indent_lvl;

  $content = 
      $indent."<table class=\"maintable\">\n".
      $indent.$tab."<tr>\n".
      $indent.$tab.$tab."<td class=\"menubase\">\n".
      $indent.$tab.$tab.$tab."<div id=\"menu\">\n".
      $indent.$tab.$tab.$tab.$tab.
      "<script src=\"../meme-suite-menu.js\" type=\"text/javascript\"></script>\n".
      $indent.$tab.$tab.$tab.$tab.
      "<script  type=\"text/javascript\">make_menu('../');</script>\n".
      $indent.$tab.$tab.$tab."</div>\n".
      $indent.$tab.$tab."</td>\n".
      $indent.$tab.$tab."<td class=\"maintablewidth\">\n".
      $indent.$tab.$tab.$tab."<div id=\"main\">\n".
      $indent.$tab.$tab.$tab.$tab."<noscript>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<h1>MEME Suite</h1>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."The MEME Suite web application requires the use of JavaScript<br />\n".
      $indent.$tab.$tab.$tab.$tab.$tab."Javascript doesn't seem to be available on your browser.\n".
      $indent.$tab.$tab.$tab.$tab."</noscript>\n".
      $indent.$tab.$tab.$tab.$tab."<form enctype=\"multipart/form-data\" method=\"POST\" action=\"$action\">\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<table>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."<col width=\"20%\">\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."<tr>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab."<td width=\"350\" valign=\"top\">\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."<img src=\"$logo\" alt=\"$alt\"><br />\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."<center>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."Version 4.9.1\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."</center>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab."</td>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab."<td align=\"left\" valign=\"bottom\">\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."<p align=\"justify\">\n".
      $description.
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab.$tab."</p>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab.$tab."</td>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."</tr>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."</table>\n";

  return($content);
} # make_submission_form_top

#
# make submission form bottom
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_submission_form_bottom {
  $logger->trace('call MemeWebUtils->make_submission_form_bottom') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content = $indent."</form>\n";

  return($content);
} # make_submission_form_bottom

# make input table
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#

# used to make sure we have a unique label each time we use the showhide feature
my $label_number;

sub make_input_table {
  $logger->trace('call MemeWebUtils->make_input_table') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($name, $left, $right, $indent_lvl, $tab, $showhide) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content = "";
  # open and close the showhide feature if needed
  my $open = "";
  my $close = "";
  # make sure we have a unique label each time we use the showhide feature
  if ($showhide) {
      # always do showhide if requested even if named undefined
      $name = "" unless defined($name);
      if (! defined($label_number)) {
	  $label_number = 0;
      } else {
	  $label_number++;
      }
      my $label = "option_".$label_number;
      $open = qq{<div id="$label\_ctrl" style="display:block;">\n<h3 class="meme">
                 <a href="javascript:show_hide('$label','$label\_ctrl')"><img src="../images/plus.png" style="border-width:0px;"/> $name</a>
                  </h3>
                </div>
                <div id="$label" style="display:none;">
                  <h3 class="meme">
                  <a href="javascript:show_hide('$label\_ctrl','$label')"><img src="../images/minus.png" style="border-width:0px;"/> $name</a>
                  </h3>};
      $close = "</div>";
  } else {
      $open = qq(<center><label class="mainheadingmedium">$name</label></center>) if $name;
      $close = "";
  }
  $content .= 
      $indent."<!-- $name args -->\n".
      "$open";
      
  $content .=
      $indent."<table>\n".
      $indent.$tab."<tr>\n".
      $indent.$tab.$tab."<td valign=\"top\" align=\"left\">\n".
      $left.
      $indent.$tab.$tab."</td>\n".
      $indent.$tab.$tab."<td valign=\"top\" align=\"right\">\n".
      $indent.$tab.$tab.$tab.
      "<!-- this is probably possible with a div, but in the absence of a css expert I'm using a table... -->\n".
      $indent.$tab.$tab.$tab."<table><tr><td align=\"left\">\n".
      $right.
      $indent.$tab.$tab.$tab."</td></tr></table>\n".
      $indent.$tab.$tab."</td>\n".
      $indent.$tab."</tr>\n".
      $indent."</table>\n";
  $content .=  # if we are using the show_hide feature, close it off
      $close;
} # make_input_table

#
# make a description field
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_description_field {
  $logger->trace('call MemeWebUtils->make_description_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($subject, $value, $indent_lvl, $tab) = @_;
  $value = "" unless defined $value;

  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;

  my $indent = $tab x $indent_lvl;

  my $content = 
      $indent."<!-- description -->\n".
      $indent."<a href=\"../help_email.html\"><b>Description</b></a> of your $subject:<br />\n".
      $indent."<input class=\"maininput\" type=\"TEXT\" size=\"40\" name=\"description\" value=\"$value\">\n";

  return($content);
} # make_description_field

#
# make address field
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_address_field {
  $logger->trace('call MemeWebUtils->make_address_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($address, $address_verify, $indent_lvl, $tab) = @_;
  $address = "" unless defined $address;
  $address_verify = "" unless defined $address_verify;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content = 
      $indent."<!-- address fields -->\n".
      $indent."Your <a href=\"../help_email.html\"><b>e-mail address:</b></a><br />\n".
      $indent."<input class=\"maininput\" type=\"TEXT\" size=\"30\" name=\"email\" value=\"$address\"><br />\n".
      $indent."Re-enter <b>e-mail address:</b><br />\n".
      $indent."<input class=\"maininput\" type=\"TEXT\" size=\"30\" name=\"email_confirm\" value=\"$address_verify\"><br />\n".
      $indent."<br />\n";

  return($content);
} # make_address_field

#
# make a motif file field
# 
# Used by
# fimo.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_motif_field {
  $logger->trace('call MemeWebUtils->make_motif_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($program, $name, $inline_name, $inline_data, $alphabet, $doc1, $doc2, $doc3, $doc4, $indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content;

  # print the input fields for the motif file unless already input
  if (! $inline_data) {
    $content = 
        $indent . "Your <a href=\"$doc1\"><b>motif</b></a> file:".
                  "<input class=\"maininput\" name=\"$name\" type=\"file\" /><br />\n".
        $indent . "<a href=\"$doc2\"><b>Sample</b></a> $doc4.<br />\n";
  } else {
    $content = 
        $indent."<p>\n".
        $indent.$tab."<H3>\n".
        $indent.$tab.$tab."$PROGRAM will search using your previously provided motif(s):\n".
        $indent.$tab.$tab."$doc3\n".
        $indent.$tab."</H3>\n";
    if ($alphabet) { # only print the alphabet if we know it as otherwise it will conflict with user selected alphabet
      $content .= $indent.$tab."<input type=\"hidden\" name=\"alphabet\" value=\"$alphabet\">\n"
    }
    $content .= 
        $indent.$tab."<input type=\"hidden\" name=\"inline_name\" value=\"$inline_name\">\n".
        $indent.$tab."<input type=\"hidden\" name=\"inline_$name\" value=\"$inline_data\">\n";

  } # no inline

  return($content);
} # make_motif_field

#
# make a list of supported databases field
# 
# Used by
# fimo.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_supported_databases_field {
  $logger->trace('call MemeWebUtils->make_supported_databases_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  #get params
  my (
    $name,          # name to give field
    $desc,          # brief db description, eg "Sequence"
    $index,         # path to the index file
    $query,         # url used to query the databases by substituting "~category~" with the category index.
    $docs,          # url for documentation on the database
    $indent_lvl,    # number of tabs to indent, zero if undefined
    $tab            # tab character, "\t" if undefined
  ) = @_;
  # set defaults for the indent
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;

  # get a list of categories from the index
  my @categories = load_categories($index);

  #
  my $content =  
      $indent."<br />$desc database to search--select <b>one</b> of the following:<br />\n".
      $indent."<br />\n".
      $indent."<!-- supported databases -->\n".
      $indent."A <a href=\"$docs\"><b>supported database</b></a>:<br />\n".
      $indent."<script type=\"text/javascript\" src=\"../selectdb.js\" defer=\"defer\"></script>\n".
      $indent."Category: <br />\n".
      $indent."<select id=\"category\" onchange=\"send_dblist_request('$query', this.value, '$name')\">\n".
      $indent.$tab."<option value=\"\" selected=\"selected\" ></option>\n";
  for (my $i = 0; $i < scalar(@categories); $i++) {
    $content .= $indent.$tab."<option value=\"$i\">" . $categories[$i] . "</option>\n";
  }
  $content .=
      $indent."</select><br />\n".
      $indent."Database:\n<br />".
      $indent."<select style=\"min-width:10em;\" id=\"$name\" name=\"$name\" disabled=\"disabled\">\n".
      $indent.$tab."<option value=\"\" selected=\"selected\"></option>\n".
      $indent."</select><br />\n".      
      $indent."<!-- attempt to fix partially cached pages causing the select boxes to be unusable -->\n".
      $indent."<script type=\"text/javascript\">\n".
      $indent.$tab."var last_db_index = 0;\n".
      $indent.$tab."function fixPage(e) {\n".
      $indent.$tab.$tab."send_dblist_request('$query', document.getElementById('category').value, '$name', last_db_index);\n".
      $indent.$tab."}\n".
      $indent.$tab."function fixPersistedPage(e) {\n".
      $indent.$tab.$tab."if (e.persisted) fixPage(e);\n".
      $indent.$tab."}\n".
      $indent.$tab."function storePageState(e) {\n".
      $indent.$tab.$tab."last_db_index = document.getElementById('$name').selectedIndex;\n".
      $indent.$tab."}\n".
      $indent.$tab."if (window.addEventListener) {\n".
      $indent.$tab.$tab."window.addEventListener('load', fixPage, false);\n".
      $indent.$tab.$tab."window.addEventListener( 'pageshow', fixPersistedPage, false );\n".
      $indent.$tab.$tab."window.addEventListener( 'pagehide', storePageState, false );\n".
      $indent.$tab."} else if (window.attachEvent) {\n".
      $indent.$tab.$tab."window.attachEvent('onload', fixPage);\n".
      $indent.$tab.$tab."window.attachEvent('onpageshow', fixPersistedPage);\n".
      $indent.$tab.$tab."window.attachEvent('onpagehide', storePageState);\n".
      $indent.$tab."} else {\n".
      $indent.$tab.$tab."window.onload = fixPage;\n".
      $indent.$tab.$tab."window.onpageshow = fixPersistedPage;\n".
      $indent.$tab.$tab."window.onpagehide = storePageState;\n".
      $indent.$tab."}\n".
      $indent."</script>\n";
  return $content;
} # make_supported_databases_field


#
# make upload FASTA field
# 
# Used by 
# fimo.pl
# glam2scan.pl
# mast.pl
# mcast.pl
#
sub make_upload_fasta_field {
  $logger->trace('call MemeWebUtils->make_upload_fasta_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($name, $maxsize) = @_;

  my $content = qq {
<!-- upload FASTA file -->
Your <a href="../doc/fasta-format.html"><b>FASTA</b></a> sequence file 
($maxsize sequence characters maximum):
<br />
<div id="upload_fasta_div">
  <input class="maininput" name="$name" type="file">
  <a onclick="clearFileInputField('upload_fasta_div')" href="javascript:noAction();"><b>Clear</b></a>
</div>
<br />
<a href="../examples/sample-kabat.seq"><b>Sample</b></a> DNA database.
  }; # end quote

  return($content);
} # make_upload_fasta_field

#
# make upload sequences field
# 
# Used by
# meme.pl
# glam2.pl
# for MEME-seq we no longer enforce maxsize but insist on
# DNA sequences: signal this by $maxsize == undef
sub make_upload_sequences_field {
  $logger->trace('call MemeWebUtils->make_upload_sequences_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($name1, $name2, $maxsize, $seq_doc, $alpha_doc,
    $format_doc,
    $filename_doc,
    $paste_doc,
    $sample_file,
    $sample_alphabet,
    $target
  ) = @_;
  # Allow target to be optional so define the default as if there was no target
  $target = "_self" unless defined $target;
  my $sizeinfo = defined($maxsize)?
  qq {The sequences may contain no more than <b>$maxsize</b>
      <a href="$alpha_doc"><b>characters</b></a>
      <br />
      total total in any of a large number of
      <a href="$format_doc"><b>formats</b></a>.} : 
  qq {The DNA sequences must be in <a href="$format_doc"><b>FASTA format</b></a>.};
  my $content = qq {
<br />
<!-- sequences fields -->
Please enter the <a href="$seq_doc">
<b>} . (defined($maxsize)?"":"DNA ") . qq{
sequences</b></a> which you believe share one or more
<br />
motifs.  $sizeinfo
<br />
<br />
Enter the <a href="$filename_doc"><b>name of a file</b></a>
containing the sequences here: 
<br />
<div id="upload_seqs_div">
  <input class="maininput" name="$name1" type="file">
  <a onclick="clearFileInputField('upload_seqs_div')" href="javascript:noAction();"><b>Clear</b></a>
</div>
<br />
<b>or</b>
<br />
the <a href="$paste_doc"><b>actual sequences</b>
</a> here (<a href="$sample_file" target="$target"><b>Sample $sample_alphabet Input Sequences</b></a>): 
<br />
<textarea name = "$name2" rows="5" cols="60"></textarea>
  }; # end quote
} # make_upload_sequences_field

#
# make radio buttons
#
# Values may be comma-separated "v1,v2"; v1 is printed.
# 
# Used by
# meme.pl
# gomo.pl
#
sub make_radio_field {
  $logger->trace('call MemeWebUtils->make_radio_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($text, $name, $checked, $values, $indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content = $indent.$text."\n";

  foreach my $value (@$values) {
    my ($v1, $v2) = split(/,/, $value);
    if ($v2 eq $checked) {
      $content .= $indent."<input class=\"mainbutton\" type=\"radio\" name=\"$name\" value=\"$v2\" checked=\"checked\" > $v1\n";
    } else {
      $content .= $indent."<input class=\"mainbutton\" type=\"radio\" name=\"$name\" value=\"$v2\" > $v1\n";
    }
  }

  return($content);
} # make_radio_field

#
# make option value list 
#
# Values may be comma-separated "v1,v2"; v2 is printed,
# and "v1,v2" is the value of the field selected.
# If there is no ",v2", then v1 is used as the printed text.
#
# Used by
# fimo.pl
# glam2scan.pl
# mast.pl
# mcast.pl
#
sub make_select_field {
  $logger->trace('call MemeWebUtils->make_select_field') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($text, $name, $selected, @values) = @_;

  my $options = "";
  foreach my $value (@values) {
    my ($v1, $v2) = split(/,/, $value); 
    $v2 = $v1 unless $v2;
    if ($value eq $selected) {
      $options .= "<option selected value='$value'> $v2\n";
    } else {
      $options .= "<option value='$value'> $v2\n";
    }
  }

  my $content = qq {
<!-- $name select field -->
$text
<select name='$name'>
$options
</select>
  }; # end quote

  return($content);
} # make_select_field

#
# make option value list 
#
# Values may be comma-separated "v1,v2"; v2 is printed,
# and "v1,v2" is the value of the field selected.
# If there is no ",v2", then v1 is used as the printed text.
#
# Used by
# fimo.pl
# glam2scan.pl
# mast.pl
# mcast.pl
#
sub make_select_field2 {
  $logger->trace('call MemeWebUtils->make_select_field2') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($text, $name, $selected, @values) = @_;

  my $options = "";
  foreach my $value (@values) {
    my ($v1, $v2) = split(/,/, $value); 
    $v2 = $v1 unless $v2;
    if ($value eq $selected) {
      $options .= "<option selected value='$v1'> $v2</option>\n";
    } else {
      $options .= "<option value='$v1'> $v2</option>\n";
    }
  }

  my $content = qq {
<!-- $name select field -->
$text
<select name='$name'>
$options
</select>
  }; # end quote

  return($content);
} # make_select_field2

#
# print submit button, reset button, email contact
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_submit_button {
  $logger->trace('call MemeWebUtils->make_submit_button') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($value, $email_contact, $indent_lvl, $tab) = @_;
  
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;

  my $indent = $tab x $indent_lvl;

  my $content = 
      $indent."<!-- submit and reset buttons -->\n".
      $indent."<center>\n".
      $indent.$tab."<input type=\"SUBMIT\" name=\"target_action\" value=\"$value\" onClick=\"return check(this.form)\">\n".
      $indent.$tab."&nbsp; &nbsp; &nbsp;\n".
      $indent.$tab."<input type=\"RESET\" value=\"Clear Input\">\n".
      $indent."</center>\n";

  return($content);
} # make_submit_button

# 
# make checkbox
#
# Note: $descr should include a </a> tag.
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# mast.pl
#
sub make_checkbox {
  $logger->trace('call MemeWebUtils->make_checkbox') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($name, $value, $descr, $checked) = @_;
  my ($check) = $checked ? " checked" : "";

  my $content = qq {
<input class="mainbutton" type="checkbox" name="$name" value="$value"$check>
$descr
  }; # end quote

  return($content);
} # make_checkbox

#
# make dna-only options
# 
# Used by
# meme.pl
# glam2.pl
# glam2scan.pl
# mast.pl
#
sub make_dna_only {
  $logger->trace('call MemeWebUtils->make_dna_only') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($options) = @_;

  my $content = qq {
<p>
<center>
<b>DNA-ONLY OPTIONS</b>
<br />
(Ignored for protein searches)
</center>
</p>
$options
  }; # end quote

  return($content);
} # make_dna_only

#
# make the complete input form
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub make_submission_form {
  $logger->trace('call MemeWebUtils->make_submission_form') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($header, $form_top, $required, $optional, $submit, $form_bottom, $tailer, $indent_lvl, $tab)  = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $indent2 = $tab x 6;
  my $content = 
      $header.
      $indent."<!-- form body -->\n".
      $indent."<body class=\"body\">\n".
      $form_top.
      $indent.$indent2."<fieldset>\n".
      $indent.$indent2.$tab."<legend>Data Submission Form</legend>\n".
      $required.
      $optional.
      $submit.
      $indent.$indent2."</fieldset>\n".
      $form_bottom.
      $tailer;

  return($content);
} # make_submission_form

#
# print the headers for a response form
# 
# Used by
# Utils
#
sub print_html_top {
  $logger->trace('call MemeWebUtils->print_html_top') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($program) = @_;

  print <<END;
<!DOCTYPE html>
<html>
  <title> $program - Verification </title>
  <body background=\"../images/bkg.jpg\">
    <hr />
END
} # print_html_top

#
# print the tailer for a response form
#
# Uses global variable $self->{NERRORS}.
# 
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub print_tailers {
  $logger->trace('call MemeWebUtils->print_tailers') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $NERRORS = $self->{NERRORS};

  # print error message summary
  if ($self->{NERRORS}) {
    my ($tobe, $booboo, $pronoun);
    if ($NERRORS == 1) {
      $tobe = "was";
      $booboo = "error";
      $pronoun = "it";
    } else {
      $tobe = "were";
      $booboo = "errors";
      $pronoun = "them";
    }
    print "</b></ol>\n";
    print "<b>There $tobe $NERRORS $booboo on the form.\n";
    print "Please correct $pronoun and try again.</b>\n";
  }

  #
  # finish the response form
  #
  print "<hr /></body></html>";
} # print_tailers

#
# Verify the job to the user via response page and email.
#
# Uses global variables $PROGRAM.
#
# Used by
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
# tomtom.pl
# spamo.pl
#
sub verify_opal_job {
  $logger->trace('call MemeWebUtils->verify_opal_job') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my (
    $result,             # return from launchJob
    $address,            # user's address
    $email_contact,          # email of support desk
    $message            # HTML text of message
  ) = @_;

  unless (eval {$result->fault}) {
    my $jobid = $result->getJobID();
    my $out_url = $result->getBaseURL();

    # make a link to the querystatus.cgi page
    my $query_url = "http://meme-suite.org//cgi-bin/querystatus.cgi";
    $query_url .= "?jobid=$jobid&service=$PROGRAM";

    # create the verification message
    my $verify = $self->make_verify_header($jobid, $query_url) . $message;

    # email the verification
    $self->send_verification($address, $email_contact, $jobid, $verify);

    # print the response form
    $self->print_html_top($PROGRAM);
    print $verify;
    print "You will also receive a confirming message at your email address: <b>$address</b>.";
  } else {
    my $code = $result->faultcode;
    my $errmsg = $result->faultstring;
    $self->whine("Your job submission resulted in a fault.<br />", "Code: $code<br />", "Message: $errmsg<br />");   
  }
  $self->print_tailers();

} # verify_opal_job

# 
# Private Method
#
# create the header for verification message
# 
# Used by
# Utils
#
sub make_verify_header {
  $logger->trace('call MemeWebUtils->make_verify_header') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($jobid, $query_url) = @_;
  my $service_url = "no";
  my $activity_url;
  if (substr $service_url, -1, 1 eq '/') {
    $activity_url = "$service_url../dashboard?command=statistics";
  }
  else {
    $activity_url = "$service_url/../dashboard?command=statistics";
  }
  my $verify = "Your job id is: <b> $jobid </b> <br />\n";
  $verify .= "You can view your job results at: <a href=\"$query_url\">$query_url</a> <br />";
  $verify .= "You can view server activity <a href=\"$activity_url\">here</a>.<br />";
  return($verify);
} # make_verify_header

#
# print an error message, bump the global error count and continue 
# Uses globals $self->{NERRORS} and $PROGRAM
#
# Used by
# Utils
# meme.pl
# fimo.pl
# glam2.pl
# glam2scan.pl
# gomo.pl
# mast.pl
# mcast.pl
#
sub whine {
  $logger->trace('call MemeWebUtils->whine') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my @msg = @_;
  $logger->debug(join(" ", @msg)) if $logger;
  #TODO make use of the multiple lines to correctly indent
  my $msg_str = join("\n", @msg);
  if ($self->{NERRORS} == 0) {
    print "Content-type: text/html\n\n";
    $self->print_html_top($PROGRAM);
    print "
      <h1>Error Report:</h1>       <hr />       <ol>
      <b>
    ";
  }
  print "<li><b>$msg_str</b>\n";
  $self->{NERRORS}++;
} # whine

#
# Private Method
#
# copy a file to standard output
#
# Used by
# Utils
#
sub copy_stdout
{
  $logger->trace('call MemeWebUtils->copy_stdout') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  open (F, "@_");
  print "<pre>";
  while (<F>) {
    print $_;
  }
  print "</pre>";
} # copy_stdout

#
# Private Method
#
# find_bad_sequence_data
# Find occurences of unknown symbols in sequence data.
#
# Returns an array of strings describing occurences of unknown symbols.
#
# Used by 
# Utils
#
sub find_bad_sequence_data {
  $logger->trace('call MemeWebUtils->find_bad_sequence_data') if $logger;
  my $self = shift;
  log_and_die("Expected Utils object") unless ref($self) eq 'MemeWebUtils';
  my $PROGRAM = $self->{PROGRAM};
  my ($bad_symbols, $fasta_data) = @_;
  my @lines = split "\n", $fasta_data;
  my $num_lines = scalar @lines;
  my (@bad_lines, $seq_name);
  for (my $i = 1; $i <= $num_lines; ++$i) {
    my $line = $lines[$i - 1];
    if ($line =~ m/>/) {
        $seq_name = $line;
        next;
    }
    if ($line =~ m/[$bad_symbols]/i) {
        push @bad_lines, "<br/>line number: $i<br/> sequence name: $seq_name<br/>sequence: $`<font color='red'>$&</font>$'\n";
    }
  }
  return @bad_lines;
} # find_bad_sequence_data

1; #modules must return true
