import sys, os, argparse

print()
print("###### Build site pages ######");
print("PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
print("# Script call: " + " ".join(sys.argv) + "\n----------");

parser = argparse.ArgumentParser(description="Gets stats from a bunch of abyss assemblies.");
parser.add_argument("--all", dest="all", help="Build all pages", action="store_true", default=False);
parser.add_argument("--index", dest="index", help="Without --all: build index.html. With --all: exlude index.html", action="store_true", default=False);
parser.add_argument("--variants", dest="variants", help="Without --all: build variants.html. With --all: exlude variants.html", action="store_true", default=False);
parser.add_argument("--iterative", dest="iterative", help="Without --all: build iterative.html. With --all: exlude iterative.html", action="store_true", default=False);
parser.add_argument("--summary", dest="summary", help="Without --all: build map_sim_summary.html. With --all: exlude map_sim_summary.html", action="store_true", default=False);
parser.add_argument("--annotations", dest="annotations", help="Without --all: build annotations.html. With --all: exlude annotations.html", action="store_true", default=False);
parser.add_argument("--people", dest="people", help="Without --all: build people.html. With --all: exlude people.html", action="store_true", default=False);
parser.add_argument("--links", dest="links", help="Without --all: build links.html. With --all: exlude links.html", action="store_true", default=False);
args = parser.parse_args();
# Input options.

#cwd = os.getcwd();
os.chdir("generators");

pages = {
    'index' : args.index,
    'variants' : args.variants,
    'iterative' : args.iterative,
    'summary' : args.summary,
    'annotations' : args.annotations,
    'people' : args.people,
    'links' : args.links
}

if args.all:
    pages = { page : False if pages[page] == True else True for page in pages };

if pages['index']:
    os.system("python index_generator.py");

if pages['variants']:
    os.system("Rscript variants_generator.r");

if pages['iterative']:
    os.system("Rscript iterative_generator.r");

if pages['summary']:
    os.system("Rscript map_sim_summary_generator.r");

if pages['annotations']:
    os.system("Rscript annotations_generator.r");

if pages['people']:
    os.system("python people_generator.py");

if pages['links']:
    os.system("python links_generator.py");
    
print("----------\nDone!");


