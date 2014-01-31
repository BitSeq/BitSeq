#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>

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

class FileHandle {
  std::FILE *fp;
  bool ispipe;
  bool isopen;
public:
  FileHandle() : fp(0), ispipe(false), isopen(false) { }

  FileHandle(std::string name)
  {
    if (name.substr(name.length()-3, 3).compare(".gz") == 0) {
      std::string tmp;

      tmp = "gunzip -c '";
      tmp += name;
      tmp += "'";
      ispipe = true;
      fp = popen(tmp.c_str(), "r");
      if (!fp)
	throw std::runtime_error("popen failed on " + tmp);
      isopen = true;
    } else {
      ispipe = false;
      fp = fopen(name.c_str(), "r");
      if (!fp)
	throw std::runtime_error("error opening file " + name);
      isopen = true;
    }
  }

  char *fgets(char * str, int size)
  {
    return std::fgets(str, size, fp);
  }

  void close()
  {
    if (isopen) {
      if (ispipe)
	pclose(fp);
      else
	fclose(fp);
      isopen = false;
    }
  }

  ~FileHandle()
  {
    close();
  }
};


class GTFEntry {
public:
  std::string seqname;
  std::string source;
  std::string feature;
  std::size_t start;
  std::size_t end;
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

    while (fp.fgets(buf, 4096) != NULL) {
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
    filepos = 0;
    return false;
  }

  FastaReader(std::string file): fp(file), buffer(""), bufferpos(0), filepos(0)
  {
    find_next_chromosome();
  }

  ~FastaReader()
  {
    fp.close();
  }

  std::string get_current_chromosome() const
  {
    return chromosome;
  }

  std::string read_string(std::size_t startpos, std::size_t endpos)
  {
    clear_buffer(startpos);
    read_to_buffer(endpos);
    return buffer.substr(startpos-bufferpos, endpos-startpos+1);
  }

private:
  void read_to_buffer(std::size_t endpos)
  {
    std::string ALPHABET="ACGTNUKSYMWRBDHV";
    char buf[BUFSIZE];
    std::string line;
    std::size_t lastchar = 0;

    while (filepos < endpos) {
      if (fp.fgets(buf, 4096) == NULL)
	throw "premature EOF";
      line = std::string(buf);
      if (line.at(0) == '>')
	throw "next chromosome reached!";
      lastchar = line.find_last_of(ALPHABET);
      buffer += line.substr(0, lastchar+1);
      filepos += lastchar+1;
    }
  }

  void clear_buffer(std::size_t startpos)
  {
    if (startpos > bufferpos) {
      if (startpos > bufferpos + buffer.length()) {
	buffer = "";
	bufferpos = filepos;
      } else {
	buffer = buffer.substr(startpos-bufferpos);
	bufferpos = startpos;
      }
    }
  }

  FileHandle fp;
  std::string chromosome;
  std::string buffer;
  std::size_t bufferpos;
  std::size_t filepos;
};


struct gtfcomp {
  bool operator() (const GTFEntry* lhs, const GTFEntry* rhs) const
  { return (lhs==0 || rhs==0 || lhs->start < rhs->start
	    || (lhs->start == rhs->start && lhs != rhs)); }
};

typedef std::set<GTFEntry *,gtfcomp> GTFSet;


class GTFStore {
public:
  GTFStore() : gtf_per_seq(), seqs() {}
  void insert(std::string chr, GTFEntry * gtf)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_seq.find(chr);
    // seq not found before
    if (it == gtf_per_seq.end()) {
      GTFSet *gtfs = new GTFSet();
      gtfs->insert(gtf);
      gtf_per_seq.insert(std::pair<std::string, GTFSet * >(chr, gtfs));
      seqs.insert(chr);
    } else {
      it->second->insert(gtf);
    }
  }

  GTFSet::iterator begin(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_seq.find(chr);
    if (it != gtf_per_seq.end())
      return it->second->begin();
    else
      throw "Seq not found";
  }

  GTFSet::iterator end(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_seq.find(chr);
    if (it != gtf_per_seq.end())
      return it->second->end();
    else
      throw "Seq not found";
  }

  std::size_t count(std::string chr)
  {
    std::map<std::string, GTFSet *>::iterator it = gtf_per_seq.find(chr);
    if (it != gtf_per_seq.end())
      return it->second->size();
    else
      return 0;
  }

  std::set<std::string>& get_seqs()
  {
    return seqs;
  }

private:
  std::map<std::string, GTFSet *> gtf_per_seq;
  std::set<std::string> seqs;
};


void write_fasta_entry_gencode(std::ostream& s, const GTFEntry *gtf,
			       std::string seq)
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
	      GTFStore& gtf_per_chromosome,
	      GTFStore& exons)
{
  GTFEntry *entry;
  FileHandle input = FileHandle(file);
  char buf[BUFSIZE];
  std::string line;

  while (input.fgets(buf, 4096) != NULL) {
    line = std::string(buf);
    if (line.at(0) != '#') {
      entry = new GTFEntry();
      entry->parse_gtf_line(line);

      if (entry->feature.compare("gene") == 0) {
	gtf.push_back(entry);
	gtf_per_chromosome.insert(entry->seqname, entry);
      } else if (entry->feature.compare("exon") == 0) {
	exons.insert(entry->extract_attribute("gene_id"), entry);
      } else {
	delete entry;
      }
    }
  }
  input.close();
}


void reconstruct_genes(GTFStore & gtf_per_chromosome, GTFStore & exons)
{
  std::set<std::string> genes = exons.get_seqs();
  std::set<std::string>::iterator gbeg, gend;
  gbeg = genes.begin();
  gend = genes.end();

  // Iterating over genes
  for (std::set<std::string>::iterator g = gbeg; g != gend ; ++g) {
    GTFSet::iterator ibeg = exons.begin(*g);
    GTFSet::iterator iend = exons.end(*g);
    std::size_t start = (*ibeg)->start;
    std::size_t end = (*ibeg)->end;
    for (GTFSet::iterator i = ibeg; i != iend ; ++i) {
      start = std::min(start, (*i)->start);
      end = std::max(end, (*i)->end);
    }
    GTFEntry *entry = new GTFEntry(*(*ibeg));
    entry->start = start;
    entry->end = end;
    std::stringstream ss;
    ss << "gene_id \"" << entry->extract_attribute("gene_id") << "\"; ";
    ss << "transcript_id \"" << entry->extract_attribute("gene_id") << "\"; ";
    ss << "gene_name \"" << entry->extract_attribute("gene_name") << "\"; ";
    ss << "gene_biotype \"" << entry->extract_attribute("gene_biotype") << "\"; ";
    ss << "transcript_name \"" << entry->extract_attribute("gene_name") << "\"; ";
    entry->attribute = ss.str();
    gtf_per_chromosome.insert(entry->seqname, entry);
  }
}


int main(int argc,char* argv[])
{
  std::vector<GTFEntry *> gtf = std::vector<GTFEntry *>();
  GTFStore gtf_per_chromosome = GTFStore();
  GTFStore exons = GTFStore();

  read_gtf(argv[1], gtf, gtf_per_chromosome, exons);
  FastaReader fasta(argv[2]);

  if (gtf.size() == 0) {
    std::cerr << "No gene entries found, reconstructing from exons."
	      << std::endl;
    reconstruct_genes(gtf_per_chromosome, exons);
  }
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

    int skipped = 0;
    for (GTFSet::iterator i=ibeg; i != iend; ++i) {
      std::string gene = (*i)->extract_attribute("gene_id");
      if (exons.count(gene) == 1) {
	//std::cerr << "[INFO] Skipping single-exon gene " << gene << std::endl;
	skipped++;
      } else {
	if ((*i)->strand == '+')
	  write_fasta_entry_gencode(std::cout,
				    *i,
				    fasta.read_string((*i)->start, (*i)->end));
	else if ((*i)->strand == '-')
	  write_fasta_entry_gencode(std::cout,
				    *i,
				    dna_reverse_complement(fasta.read_string((*i)->start, (*i)->end)));
	else
	  std::cerr << "Unknown strand" << std::endl;
      }
    }
    if (skipped > 0) {
      std::cerr << "[INFO] " << skipped << " single exon genes skipped."
		<< std::endl;
      skipped = 0;
    }
  } while (fasta.find_next_chromosome());

  return 0;
}
