#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>

const int BUFSIZE=4096;

std::string dna_reverse_complement(std::string orig)
{
  std::string res = std::string(orig.length(), ' ');
  std::string::iterator i=res.begin();
  std::string::reverse_iterator j=orig.rbegin();
  while (j != orig.rend()) {
    switch (*j) {
      case 'A':
	*i = 'T'; break;
      case 'C':
	*i = 'G'; break;
      case 'G':
	*i = 'C'; break;
      case 'T':
	*i = 'A'; break;
      default:
	*i = *j;
      }
    ++i;
    ++j;
  }

  return res;
}

class GTFEntry {
public:
  std::string seqname;
  std::string source;
  std::string feature;
  int start;
  int end;
  double score;
  char strand;
  char frame;
  std::string attribute;

  int parse_gtf_line(const std::string& line)
  {
    std::size_t startpos = 0;
    std::size_t endpos = line.find('\t');
    seqname = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    source = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    feature = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    start = std::atoi(line.substr(startpos, endpos-startpos).c_str());
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    end = std::atoi(line.substr(startpos, endpos-startpos).c_str());
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    std::string tmp = line.substr(startpos, endpos-startpos);
    if (tmp.compare(std::string(1, '.')) != 0)
      score = std::atoi(tmp.c_str());
    else
      score = 0.0;
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    strand = line.at(startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    frame = line.at(startpos);
    startpos = endpos + 1;

    endpos = line.find('\t', startpos);
    attribute = line.substr(startpos, endpos-startpos);
    startpos = endpos + 1;

    return 0;
  }

  std::string extract_attribute(std::string name) const
  {
    std::size_t startpos = 0;
    std::size_t endpos = attribute.find(' ');
    while (endpos < std::string::npos) {
      if (attribute.compare(startpos, endpos, name) == 0) {
	std::size_t terminatorpos = attribute.find(';', endpos);
	return attribute.substr(endpos+1, terminatorpos-endpos-1);
      }
    }
    return std::string("");
  }

  void pretty_print(std::ostream& s) const
  {
    s << "seqname: " << seqname << std::endl;
    s << "source: " << source << std::endl;
    s << "feature: " << feature << std::endl;
    s << "start: " << start << std::endl;
    s << "end: " << end << std::endl;
    s << "score: " << score << std::endl;
    s << "strand: " << strand << std::endl;
    s << "frame: " << frame << std::endl;
    s << "attribute: " << attribute << std::endl;
  }
};


class FastaReader {
public:
  void find_next_chromosome()
  {
    char buf[BUFSIZE];
    std::string line;

    while (std::fgets(buf, 4096, fp) != NULL) {
      line = std::string(buf);
      if (line.at(0) == '>') {
	std::size_t endpos = line.find(' ', 1);
	chromosome = line.substr(1, endpos-1);
	filepos = 1;
	buffer = "";
	bufferpos = 1;
	return;
      }
    }
    chromosome = "";
    filepos = -1;
  }

  FastaReader(std::string file): buffer(""), bufferpos(0), filepos(0)
  {
    std::string tmp;

    tmp = "gunzip -c '";
    tmp += file;
    tmp += "'";
    fp = popen(tmp.c_str(), "r");
    find_next_chromosome();
  }

  ~FastaReader()
  {
    fclose(fp);
  }

  std::string get_current_chromosome() const
  {
    return chromosome;
  }

  std::string read_string(int startpos, int endpos)
  {
    clear_buffer(startpos);
    read_to_buffer(endpos);
    return buffer.substr(startpos-bufferpos, endpos-startpos+1);
  }

private:
  void read_to_buffer(int endpos)
  {
    std::string ALPHABET="ACGTNUKSYMWRBDHV";
    char buf[BUFSIZE];
    std::string line;

    while (std::fgets(buf, 4096, fp) != NULL && filepos < endpos) {
      line = std::string(buf);
      std::size_t endpos = line.find_last_of(ALPHABET);
      buffer += line.substr(0, endpos+1);
      filepos += endpos;
    }
  }

  void clear_buffer(int startpos)
  {
    if (startpos > bufferpos) {
      if (startpos-bufferpos > buffer.length()) {
	buffer = "";
	bufferpos = filepos;
      } else {
	buffer = buffer.substr(startpos-bufferpos);
	bufferpos = startpos;
      }
    }
  }

  std::FILE *fp;
  std::string chromosome;
  std::string buffer;
  int bufferpos;
  int filepos;
};


void read_gtf(std::string file, std::vector<GTFEntry *>& gtf,
	      std::multimap<std::string, GTFEntry *>& gtf_per_chromosome,
	      std::set<std::string>& chromosomes)
{
  GTFEntry *entry;
  std::FILE *input;
  char buf[BUFSIZE];
  std::string tmp;
  std::string line;

  tmp = "gunzip -c '";
  tmp += file;
  tmp += "'";
  input = popen(tmp.c_str(), "r");
  while (std::fgets(buf, 4096, input) != NULL) {
    line = std::string(buf);
    if (line.at(0) != '#') {
      entry = new GTFEntry();
      entry->parse_gtf_line(line);

      if (entry->feature.compare("gene") == 0) {
	gtf.push_back(entry);
	gtf_per_chromosome.insert(std::pair<std::string, GTFEntry *>(entry->seqname, entry));
	chromosomes.insert(entry->seqname);
      } else {
	delete entry;
      }
    }
  }
  fclose(input);
}


int main(int argc,char* argv[])
{
  std::vector<GTFEntry *> gtf = std::vector<GTFEntry *>();
  std::multimap<std::string, GTFEntry *> gtf_per_chromosome =
    std::multimap<std::string, GTFEntry *>();
  std::set<std::string> chromosomes = std::set<std::string>();

  read_gtf(argv[1], gtf, gtf_per_chromosome, chromosomes);
  FastaReader fasta(argv[2]);

  gtf[0]->pretty_print(std::cout);
  std::cout << "gene_id: " << gtf[0]->extract_attribute("gene_id")
	    << std::endl;
  std::cout << gtf.size() << " genes found." << std::endl;
  for (std::set<std::string>::const_iterator i=chromosomes.begin();
       i != chromosomes.end(); i++) {
    std::cout << *i << ": " << gtf_per_chromosome.count(*i) 
	      << " genes." << std::endl;
  }

  std::cout << "Fasta chromosome: " << fasta.get_current_chromosome()
	    << std::endl;

  std::pair<std::multimap<std::string, GTFEntry *>::iterator, std::multimap<std::string, GTFEntry *>::iterator> range;
  range = gtf_per_chromosome.equal_range(fasta.get_current_chromosome());

  for (std::multimap<std::string, GTFEntry *>::iterator i=range.first;
       i != range.second; i++) {
    std::cout << i->second->start << " " << i->second->end << std::endl;
    if (i->second->strand == '+')
      std::cout << fasta.read_string(i->second->start, i->second->end) << std::endl;
    else if (i->second->strand == '-')
      std::cout << dna_reverse_complement(fasta.read_string(i->second->start, i->second->end))
	      << std::endl;
    else
      std::cout << "Unknown strand" << std::endl;
  }

  return 0;
}
