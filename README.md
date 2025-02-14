```
Παπαλέξης Γεράσιμος 1115202000260
Παπαλέξης Σπυρίδων 1115201900148

Github link: https://github.com/SpyrosPapalexis/Project01

Το πρόγραμμα που υλοποιήσαμε διαβάζει τα δεδομένα εισόδου από ένα json αρχείο το οποίο διαβάζεται είτε ως όρισμα στην κλήση του, είτε χειροκίνητα αν δεν δωθεί από τον χρήστη.
Στο αρχείο json καταγράφεται η μεθοδολογία τριγωνοποίησης, η οποία μπορεί να είναι είτε local, είτε sa, είτε ant ή brute-force (υπολοίηση πρώτης εργασίας) στην περίπτωση που η μεταβλητή delaunay είναι ψευδής.
    1. τοπική αναζήτηση: ελέγχει κάθε μέθοδο σε κάθε κύκλο και επιλέγει αυτή που βελτιστοποιεί την τριγωνοποίηση (έχει στο πλήθος λιγότερα αμβλυγώνια).

    2. προσομοιωμένη ανόπτηση: επιλέγει μέθοδο εισαγωγής τυχαία για κάθε steiner και υπολογίζει την ενέργειά της. Η μέθοδος γίνεται αποδεκτή αν έχει λιγότερη ενέργεια ή με το κριτήριο Metropolis.

    3. αποικία μυρμηγκιών: το κάθε μυρμήγκι επιλέγει μια μέθοδο χρησιμοποιώντας ως βάρος τις ευρετικές και τις φερορμόνες. Στο τέλος κάθε κύκλου οι φερορμόνες ενημερώνονται ανάλογα με το ποιές μεθόδοι βοήθησαν στην τριγωνοποίηση.

    4. brute-force: ελέγχει κάθε μέθοδο εισαγωγής για L steiner points και επιλέγει αυτήν με το ελάχιστο πλήθος αμβλυγώνιων τριγώνων.

Οι μέθοδοι τριγωνοποίησης βρίσκουν και επιλέγουν το πρώτο αμβλυγώνιο τρίγωνο το οποίο είναι εντός του boundary. Όταν η εισαγωγή συμπίπτει σε ήδη υπαρκτό steiner ή σε constrained edge, επιλέγει το επόμενο τρίγωνο. Στην περίπτωση όπου το τρίγωνο είναι nullptr, σημαίνει ότι δεν υπάρχει άλλο αμβλυγώνιο τρίγωνο και κάθε άρα σταματάει την εισαγωγή.
Οι μέθοδοι εισαγωγής steiner είναι οι εξής:
    1. steiner point στο μέσο της μεγαλύτερης πλευράς ενός αμβλυγώνιου τριγώνου.
    2. steiner point στο περίκεντρο ενός αμβλυγώνιου τριγώνου. (ή στο βαρύκεντρο, σε περίπτωση που το περίκεντρο βρίσκεται εκτός του boundary)
    3. steiner point στο σημείο τομής του ύψους της προβολής μιας αμβλίας γωνίας με την προβολή αυτή.
    4. steiner point στο κεντροειδές σημείο ενός αμβλυγώνιου τριγώνου.

Το πρόγραμμα στην αρχή εκτυπώνει το αρχικό πλήθος αμβλυγώνιων τριγώνων, καθώς και την αρχική οπτικοποίηση. Ενώ στην συνέχεια εκτυπώνει το τελικό πλήθος αμβλυγώνιων τριγώνων, εισαχθέντων steiner, καθώς και την τελική οπτικοποίηση.
Το πρόγραμμα αποτελείται από ένα αρχείο με όνομα "cgalexec.cpp".
Η έξοδος ονομάζεται "output.json".

Το πρόγραμμα μεταγλωτίζεται με την κλήση:
cmake -DCGAL_DIR=/usr/lib/CGAL
make