#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <exception>

const int BUFSIZE=4096;
const unsigned int LINEWIDTH=70;

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

inline std::string dequote(std::string orig) {
  if (orig.at(0) == '"' && orig.at(orig.length()-1) == '"')
    return orig.substr(1, orig.length()-2);
  else
    return orig;
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
    std::size_t terminatorpos;
    while (endpos < std::string::npos) {
      terminatorpos = attribute.find(';', endpos);
      if (attribute.compare(startpos, endpos-startpos, name) == 0) {
	return dequote(attribute.substr(endpos+1, terminatorpos-endpos-1));
      }
      startpos = terminatorpos + 2;
      endpos = attribute.find(' ', startpos);
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
  bool find_next_chromosome()
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
	return true;
      }
    }
    chromosome = "";
    filepos = -1;
    return false;
  }

  FastaReader(std::string file): buffer(""), bufferpos(0), filepos(0)
  {
    std::string tmp;

    tmp = "gunzip -c '";
    tmp += file;
    tmp += "'";
    fp = popen(tmp.c_str(), "r");
    if (!fp)
      throw std::runtime_error("popen failed on " + tmp);
    find_next_chromosome();
  }

  ~FastaReader()
  {
    pclose(fp);
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
    std::size_t lastchar = 0;

    while (filepos < endpos) {
      if (std::fgets(buf, 4096, fp) == NULL)
	throw "premature EOF";
      line = std::string(buf);
      if (line.at(0) == '>')
	throw "next chromosome reached!";
      lastchar = line.find_last_of(ALPHABET);
      buffer += line.substr(0, lastchar+1);
      filepos += lastchar+1;
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


struct gtfcomp {
  bool operator() (const GTFEntry* lhs, const GTFEntry* rhs) const
  { return (lhs==0 || rhs==0 || lhs->start < rhs->start
	    || (lhs->start == rhs->start && lhs != rhs)); }
};

typedef std::set<GTFEntry *,gtfcomp> GTFSet;


class GTFStore {
public:
  GTFStore() : gtf_per_chromosome(), chromosomes() {}
  void insert(std::string chr, GTFEntry * gtf)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_chromosome.find(chr);
    // chromosome not found before
    if (it == gtf_per_chromosome.end()) {
      GTFSet *gtfs = new GTFSet();
      gtfs->insert(gtf);
      gtf_per_chromosome.insert(std::pair<std::string, GTFSet * >(chr, gtfs));
      chromosomes.insert(chr);
    } else {
      it->second->insert(gtf);
    }
  }

  GTFSet::iterator begin(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_chromosome.find(chr);
    if (it != gtf_per_chromosome.end())
      return it->second->begin();
    else
      throw "Chromosome not found";
  }

  GTFSet::iterator end(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_chromosome.find(chr);
    if (it != gtf_per_chromosome.end())
      return it->second->end();
    else
      throw "Chromosome not found";
  }

  std::size_t count(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_chromosome.find(chr);
    if (it != gtf_per_chromosome.end())
      return it->second->size();
    else
      return 0;
  }

  std::set<std::string>& get_chromosomes()
  {
    return chromosomes;
  }

private:
  std::map<std::string, GTFSet *> gtf_per_chromosome;
  std::set<std::string> chromosomes;
};


void write_fasta_entry(std::ostream& s, const GTFEntry *gtf, std::string seq)
{
  std::string tmp;

  s << ">";
  s << gtf->extract_attribute("transcript_id") << "|";
  s << gtf->extract_attribute("gene_id") << "|";

  tmp = gtf->extract_attribute("havana_gene");
  if (tmp.empty())
    s << "-|";
  else
    s << tmp << "|";

  tmp = gtf->extract_attribute("havana_transcript");
  if (tmp.empty())
    s << "-|";
  else
    s << tmp << "|";
  
  s << gtf->extract_attribute("transcript_name") << "|";
  s << gtf->extract_attribute("gene_name") << "|";
  s << seq.length() << "|" << std::endl;

  int start = 0;
  while (seq.length() > start + LINEWIDTH) {
    s << seq.substr(start, LINEWIDTH) << std::endl;
    start += LINEWIDTH;
  }
  s << seq.substr(start) << std::endl;
}


void read_gtf(std::string file, std::vector<GTFEntry *>& gtf,
	      GTFStore& gtf_per_chromosome)
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
	gtf_per_chromosome.insert(entry->seqname, entry);
      } else {
	delete entry;
      }
    }
  }
  pclose(input);
}


int main(int argc,char* argv[])
{
  std::vector<GTFEntry *> gtf = std::vector<GTFEntry *>();
  GTFStore gtf_per_chromosome = GTFStore();

  read_gtf(argv[1], gtf, gtf_per_chromosome);
  FastaReader fasta(argv[2]);

  std::cerr << gtf.size() << " genes found." << std::endl;

  GTFSet::iterator ibeg;
  GTFSet::iterator iend;
  do {
    std::string chr = fasta.get_current_chromosome();
    
    try {
      ibeg = gtf_per_chromosome.begin(chr);
      iend = gtf_per_chromosome.end(chr);
    } catch (...) {
      std::cerr << "No gtf entries for chromosome " << chr << std::endl;
      continue;
    }
    std::cerr << "Starting chromosome " << chr << " with " << gtf_per_chromosome.count(chr) << " genes." << std::endl;

    for (GTFSet::iterator i=ibeg; i != iend; ++i) {
      if ((*i)->strand == '+')
	write_fasta_entry(std::cout,
			  *i,
			  fasta.read_string((*i)->start, (*i)->end));
      else if ((*i)->strand == '-')
	write_fasta_entry(std::cout,
			  *i,
			  dna_reverse_complement(fasta.read_string((*i)->start, (*i)->end)));
      else
	std::cerr << "Unknown strand" << std::endl;
    }
  } while (fasta.find_next_chromosome());

  return 0;
}
