from collections import Counter
import pandas as pd
from os.path import join
import numpy as np


# assign all genera to taxonomic families
families = {
    'Octopodidae': ['Octopus'],
    'Scophthalmidae': ['Scophthalmus'],
    'Acipenseridae': ['Acipenser'],
    'Sciaenidae': ['Argyrosomus'],
    'Parastacidae': ['Cherax'],
    'Clariidae': ['Clarias'],
    'Cyprinidae': ['Ctenopharyngodon', 'Cyprinus', 'Pethia', 'Puntius', 'Labeo', 'Hypophthalmichthys'],
    'Sparidae': ['Dentex', 'Diplodus', 'Pagrus'],
    'Unknown': ['Dicentrachus'],
    'Serranidae': ['Epinephelus'],
    'Gadidae': ['Gadus'],
    'Pleuronectidae': ['Hippoglossus'],
    'Nephropidae': ['Homarus'],
    'Latidae': ['Lates', 'Lota'],
    'Penaeidae': ['Litopenaeus'],
    'Lutjanidae': ['Lutjanus'],
    'Moronidae': ['Morone'],
    'Mugilidae': ['Mugil'],
    'Salmonidae': ['Oncorhynchus', 'Salvelinus', 'Thymallus', 'Salmo'],
    'Cichlidae': ['Oreochromis', 'Pterophyllum', 'Neolamprologus'],
    'Pangasiidae': ['Pangasianodon'],
    'Penaeidae': ['Penaeus'],
    'Percidae': ['Perca', 'Sander'],
    'Polyprionidae': ['Polyprion'],
    'Rachycentridae': ['Rachycentron'],
    'Sepiidae': ['Sepia'],
    'Carangidae': ['Seriola', 'Trachinotus'],
    'Soleidae': ['Solea'],
    'Sparidae': ['Sparus'],
    'Scombridae': ['Thunnus'],
    'Osphronemidae': ['Betta'],
    'Cyprinidae': ['Carassius', 'Danio', 'Megalobrama', 'Mylopharyngodon', 'Cirrhinus', 'Barbonymus'],
    'Poeciliidae': ['Poecilia'],
    'Characidae': ['Hyphessobrycon'],
    'Callichthyidae': ['Corydoras'],
    'Botiidae': ['Chromobotia'],
    'Osphronemidae': ['Osphronemus'],
    'Channidae': ['Channa', 'Chanos'],
    'Varunidae': ['Eriocheir'],
    'Cambaridae': ['Procambarus'],
    'Acipenseridae': ['Huso'],
    'Mugilidae': ['Chelon'],
    'Centrarchidae': ['Micropterus'],
    'Sinipercidae': ['Siniperca'],
    'Lateolabracidae': ['Lateolabrax'],
    'Sciaenidae': ['Larimichthys', 'Sciaenops'],
    'Bagridae': ['Tachysurus'],
    'Palaemonidae': ['Macrobrachium'],
    'Ictaluridae': ['Ictalurus'],
    'Siluridae': ['Silurus'],
    'Cobitidae': ['Misgurnus'],
    'Synbranchidae': ['Monopterus'],
    'Coryphaenidae': ['Coryphaena']
}

# assign all families to taxonomic classes
classes = {
    'Cephalopoda': ['Octopodidae', 'Sepiidae'],
    'Actinopterygii': ['Scophthalmidae', 'Acipenseridae', 'Sciaenidae',
                       'Clariidae', 'Cyprinidae', 'Sparidae', 'Serranidae',
                       'Gadidae', 'Pleuronectidae', 'Latidae', 'Lutjanidae',
                       'Moronidae', 'Mugilidae', 'Salmonidae', 'Cichlidae',
                       'Pangasiidae', 'Percidae', 'Polyprionidae', 'Rachycentridae',
                       'Carangidae', 'Soleidae', 'Scombridae', 'Osphronemidae',
                       'Poeciliidae', 'Characidae', 'Callichthyidae',
                       'Botiidae', 'Channidae', 'Centrarchidae', 'Sinipercidae',
                       'Lateolabracidae', 'Bagridae', 'Ictaluridae', 'Siluridae',
                       'Cobitidae', 'Synbranchidae', 'Coryphaenidae'],
    'Malacostraca': ['Parastacidae', 'Nephropidae', 'Penaeidae', 'Varunidae',
                     'Cambaridae', 'Palaemonidae'],
}

wd = join('..', '..', 'data')
df = pd.read_csv(join(wd, 'feb_farm_data_s2.csv'))

species = [*df['Name'].unique()]
species.remove(np.nan)
genera_repeated = [x.split(' ')[0] for x in species]
genera = [*set(genera_repeated)]

genus_counts = Counter(genera_repeated)
print(genus_counts.most_common()[:3])  # show three most common genera
genus_counts = dict(genus_counts)

family_counts = {}
for family_i, genera_i in families.items():
        count = sum([genus_counts[genus] for genus in genera_i])
        family_counts[family_i] = count

class_counts = {}
for class_i, families_i in classes.items():
        count = sum([family_counts[family] for family in families_i])
        class_counts[class_i] = count


# create dictionaries assigning genera, then species to classes
classes_genera = {
    'Cephalopoda': [],
    'Actinopterygii': [],
    'Malacostraca': []
}
classes_species = {}
for class_i, families_i in classes.items():
    for family in families_i:
        classes_genera[class_i] += families[family]

for class_i, genera_i in classes_genera.items():
    classes_species[class_i] = []

    for genus in genera_i:
        classes_species[class_i] += [*filter(lambda x: genus in x, species)]

for class_i, species_i in classes_species.items():
    print(f"{class_i}:, {len(species_i)}")
