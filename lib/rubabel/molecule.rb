require 'open-uri'
require 'rexml/document'
require 'openbabel'
require 'rubabel'
require 'rubabel/atom'
require 'rubabel/bond'
require 'stringio'
require 'mini_magick'
require 'set' # Has one downside: only stores unique entries, which might be a problem later

class OpenBabel::OBMol
  def upcast
    Rubabel::Molecule.new(self)
  end
end

class OpenBabel::OBConversion
  OUT_OPTS_SHORT = {
    no_element_coloring: :u,
    no_internal_specified_color: :U,
    black_bkg: :b,
    no_terminal_c: :C,
    draw_all_carbon: :a,
    no_molecule_name: :d,
    embed_mol_as_cml: :e,
    scale_to_bondlength_pixels: :p,
    size: :P, # single number of pixels (e.g., 300, but will take dims for png: ('300x250'))
    add_atom_index: :i,
    no_javascript: :j,
    wedge_hash_bonds: :w,
    no_xml_declaration: :x,
    aliases: :A,
  }

  # adds the opts to the type, where type is :gen, :out, or :in.
  # takes opts as either a hash or a list:
  #
  #     # hash
  #     obconv.add_opts!(:out, d: true, u: true, p: 10 ) 
  #
  #     # list
  #     obconv.add_opts!(:out, :d, :u, [:p, 10] )
  #
  # returns self
  def add_opts!(type=:gen, *opts)
    opt_type = OpenBabel::OBConversion.const_get(type.to_s.upcase.<<("OPTIONS"))
    hash = 
      if opts.first.is_a?(Hash)
        opts.first
      else
        opts.inject({}) do |hsh,v| 
          if v.is_a?(Array)
            hsh[v[0]] = v[1]
          else
            hsh[v] = true
          end
          hsh
        end
      end
    hash.each do |k,v|
      next if v == false
      args = [ (OUT_OPTS_SHORT[k] || k).to_s, opt_type ]
      if v && (v != true)
        args << v
      end
      self.add_option(*args)
    end
    self
  end
end

class OpenBabelUnableToSetupForceFieldError < RuntimeError
end

module Rubabel
  class Molecule
    include Enumerable

    DEFAULT_FINGERPRINT = "FP2"
    DEFAULT_OUT_TYPE = :can
    DEFAULT_IN_TYPE = :smi
    DEFAULT_ID_TYPE = :pubchem
    ID_TYPES = {
      inchikey: "InChI key",
      lmid: "LipidMaps ID"
    }
    ID_TYPE_KEYS = ID_TYPES.keys
    ADDUCTS_LEGEND = {
      lithium: "[Li+]",
      sodium: "[Na+]",
      ammonium: "[NH4+]",
      potassium: "[K+]",
      hydrogen: "[H+]",
      #name: "[ion+]",
      hydride: "[H-]",
      hydroxide: "[OH-]",
      chlorine: "[Cl-]",
      bromine: "[Br-]" #,
      #name2: "[ion-]"
    }
    ADDUCTS = {}
    
    @adducts = []
    # the OpenBabel::OBmol object is :ob
    attr_accessor :ob
    attr_accessor :adducts


    class << self

      def tanimoto(mol1, mol2, type=DEFAULT_FINGERPRINT)
        OpenBabel::OBFingerprint.tanimoto(mol1.ob_fingerprint(type), mol2.ob_fingerprint(type))
      end

      def from_file(file, type=nil)
        (obmol, obconv, not_at_end) = Rubabel.read_first_obmol(file, type).first
        Rubabel::Molecule.new(obmol)
      end

      def load_archive
        @archive = YAML.load_file(LMID_ARCHIVE)
      end
      def from_archive(lmid)
        if Rubabel::ARCHIVE # archive is turned on
          load_archive unless @archive
          if @archive[lmid]
            # Lookup from the archive based on the input string (AN LMID)
            string = @archive[lmid][:structure]
            unless lmid = @archive[lmid][:lmid]
              lmid = lmid[/LM[A-Z]{2}\d{8,10}/]
            end
            type = :sdf
            obmol = OpenBabel::OBMol.new
            obconv = OpenBabel::OBConversion.new
            obconv.set_in_format(type.to_s) || raise(ArgumentError, "invalid format #{type}")
            obconv.read_string(obmol, string) || raise(ArgumentError, "invalid string" )
            resp = self.new(obmol)
          else
            unless lmid.size < 12
              lmid = lmid[/LM[A-Z]{2}\d{8,10}/]
            end
            url = "http://www.lipidmaps.org/data/LMSDRecord.php?OutputType=SDF&Mode=File&LMID=" + lmid
            doc_string = retrieve_info_from_url(url)
            resp = from_string(doc_string, :sdf)
          end
        else
          url = "http://www.lipidmaps.org/data/LMSDRecord.php?OutputType=SDF&Mode=File&LMID=" + lmid
          doc_string = retrieve_info_from_url(url)
          resp = from_string(doc_string, :sdf)
        end
      end




      def from_string(string, type=DEFAULT_IN_TYPE, adducts=nil)
        if type == :inchi
          string.prepend("InChI=") unless string[/^InChI=/]
        end
        obmol = OpenBabel::OBMol.new
        obconv = OpenBabel::OBConversion.new
        obconv.set_in_format(type.to_s) || raise(ArgumentError, "invalid format #{type}")
        obconv.read_string(obmol, string) || raise(ArgumentError, "invalid string" )
        ob = self.new(obmol)
        if adducts
          ob.adducts = Array(adducts)
          ob.adducts = ob.adducts.map do |m| 
            from_string(m) if m.class != Rubabel::Molecule
          end
        else
          ob.find_internal_adducts(ob)
        end
        ob
      end

      def retrieve_info_from_url(url)
        begin
          info = open(url) {|io| io.read }
        rescue => e
          puts "Some kind of internet connectivity error. Check your connection!"
          raise e
        end
        info
      end

      # requires an internet connection
      def from_id(id, type=DEFAULT_ID_TYPE)
          case type
          when :inchikey
            url = "http://www.chemspider.com/InChI.asmx/InChIKeyToInChI?inchi_key=" + URI::encode(id)
            doc_string = retrieve_info_from_url(url)
            doc = REXML::Document.new( doc_string )
            inchi_string = doc.root.children.first.to_s
            raise(ArgumentError, "did not retrieve a valid inchi string") unless inchi_string[/^InChI=/]
            from_string(inchi_string, :inchi)
          when :lmid # lipidmaps id
            url = "http://www.lipidmaps.org/data/LMSDRecord.php?OutputType=SDF&Mode=File&LMID=" + id
            doc_string = retrieve_info_from_url(url)
            from_string(doc_string, :sdf)
          end
      end
      alias_method :from_key, :from_id

      def from_atoms_and_bonds(atoms=[], bonds=[])
        obj = self.new( OpenBabel::OBMol.new )
        atoms.each {|atom| obj.add_atom(atom) }
        bonds.each {|bond| obj.add_bond(bond) }
        obj
      end
    end


    #Code which can be used to find addcuts from molecule containing them
    def find_internal_adducts(mol = self)
      mol.adducts ||= []
      splits = mol.ob.separate.map(&:upcast).flatten.each {|m| mol.adducts << m if ADDUCTS.values.include?(m) }
      mol.adducts = mol.adducts.uniq(&:csmiles)
      # TODO add a filter to remove the captured adducts from the molecule
      #??? mol = mol.ob.separate.map(&:upcast).delete_if {|m| ADDUCTS.values.include?(m) }
      mol.adducts = mol.adducts.map {|m| m.class == Rubabel::Molecule ? m : from_string(m)  }
    end
    def adducts
      @adducts ||= []
      @adducts = @adducts.map {|m| m.class == Rubabel::Molecule ? m : Rubabel::Molecule.from_string(m)  }
      @adducts = @adducts.uniq(&:csmiles)
    end
    def adduct?
      self.find_internal_adducts unless @adducts
      @adducts.size > 0
    end
    alias :adducts? :adduct?
    alias :adduct :adducts



    # arg may be a Fixnum, a Symbol (Elemental symbol that is a Symbol), or a
    # Rubabel::Atom.  Returns the newly associated/created atom.
    def associate_atom!(arg)
      if arg.is_a?(Rubabel::Atom)
        @ob.add_atom(arg.ob)
        arg
      else
        (num, is_aromatic) = 
          if arg.is_a?(Symbol)
            [Rubabel::ELEMENT_TO_NUM[arg], (arg.to_s.capitalize != arg.to_s)]
          else
            [arg, false]
          end

        new_obatom = @ob.new_atom
        new_obatom.set_atomic_num(num)
        new_obatom.set_aromatic if is_aromatic
        Rubabel::Atom.new(new_obatom)
      end
    end


    # returns the atom passed in or that was created.  arg is a pre-existing
    # atom, an atomic number or an element symbol (e.g. :c).  default is to
    # add carbon.
    def add_atom!(arg=6, bond_order=1, attach_to=nil)
      attach_to ||= atoms.last
      atom = associate_atom!(arg)
      add_bond!(attach_to, atom, bond_order) if attach_to
      atom
    end

    alias_method :<<, :add_atom!

    # retrieves the atom by index (accepts everything an array would)
    def [](*args)
      atoms[*args]
    end

    def delete_atom(atom)
      @ob.delete_atom(atom.ob, false)
    end

    # attributes
    def title() @ob.get_title end
    def title=(val) @ob.set_title(val) end

    def charge() @ob.get_total_charge end
    def charge=(v) @ob.set_total_charge(v) end

    def spin() @ob.get_total_spin_multiplicity end

    def mol_wt() @ob.get_mol_wt end
    alias_method :avg_mass, :mol_wt

    def exact_mass() @ob.get_exact_mass end

    # returns the exact_mass corrected for charge gain/loss
    def mass
      @ob.get_exact_mass - (@ob.get_total_charge * Rubabel::MASS_E)
    end

    def mass_with_adduct
      [mass].product(@adducts.map(&:mass)).map{|a| a.inject(:+) }
    end
    def masses_with_possible_adducts
      [mass].product([0,@adducts.map(&:mass)].flatten).map {|a| a.inject(:+) }
    end

    # returns a string representation of the molecular formula.  Not sensitive
    # to add_h!
    def formula() @ob.get_formula end

    def initialize(obmol=nil)
      @ob = obmol.nil? ? OpenBabel::OBMol.new  : obmol
    end

    # returns a list of atom indices matching the patterns (corresponds to the
    # OBSmartsPattern::GetUMapList() method if uniq==true and GetMapList
    # method if uniq==false).  Note that the original GetUMapList returns atom
    # *numbers* (i.e., the index + 1).  This method returns the zero indexed
    # indices.
    def smarts_indices(smarts_or_string, uniq=true)
      mthd = uniq ? :get_umap_list : :get_map_list
      pattern = smarts_or_string.is_a?(Rubabel::Smarts) ? smarts_or_string : Rubabel::Smarts.new(smarts_or_string)
      pattern.ob.match(@ob)
      pattern.ob.send(mthd).map {|atm_indices| atm_indices.map {|i| i - 1 } }
    end

    # yields atom arrays matching the pattern.  returns an enumerator if no
    # block is given
    def each_match(smarts_or_string, uniq=true, &block)
      block or return enum_for(__method__, smarts_or_string, uniq)
      _atoms = self.atoms
      smarts_indices(smarts_or_string, uniq).each do |ar|
        block.call(_atoms.values_at(*ar))
      end
    end

    # returns an array of matching atom sets.  Consider using each_match.
    def matches(smarts_or_string, uniq=true)
      each_match(smarts_or_string, uniq).map.to_a
    end

    def matches?(smarts_or_string)
      # TODO: probably a more efficient way to do this using API
      smarts_indices(smarts_or_string).size > 0
    end

    # returns an array of OpenBabel::OBRing objects.
    def ob_sssr
      @ob.get_sssr.to_a
    end

    #def conformers
    # Currently returns an object of type
    # SWIG::TYPE_p_std__vectorT_double_p_std__allocatorT_double_p_t_t
    #vec = @ob.get_conformers
    #end

    # are there hydrogens added yet
    def hydrogens_added?
      @ob.has_hydrogens_added
    end
    alias_method :h_added?, :hydrogens_added?

    # ensures that hydrogens are added before an operation, but returns the
    # molecule to the original hydrogen or no hydrogen state when finished.
    # returns whatever was returned by the block.
    def do_with_hydrogens(ph=nil, polaronly=false, &block)
      hydr_added = @ob.has_hydrogens_added
      add_h!(ph, polaronly) unless hydr_added
      reply = block.call
      @ob.delete_hydrogens unless hydr_added
      reply
    end

    def do_without_hydrogens(ph=nil, polaronly=false, &block)
      hydr_added = @ob.has_hydrogens_added
      @ob.delete_hydrogens if hydr_added
      reply = block.call
      add_h!(ph, polaronly) if hydr_added
      reply
    end

    # returns self.  Corrects for ph if ph is not nil.  NOTE: the reversal of
    # arguments from the OpenBabel api.
    def add_h!(ph=nil, polaronly=false)
      if ph.nil?
        @ob.add_hydrogens(polaronly)
      else
        @ob.add_hydrogens(polaronly, true, ph)
      end
      self
    end

    # only adds polar hydrogens. returns self
    def add_polar_h!
      @ob.add_polar_hydrogens
      self
    end

    # returns self.  If ph is nil, then #neutral! is called
    def correct_for_ph!(ph=7.4)
      ph.nil? ? neutral! : @ob.correct_for_ph(ph)
      self
    end

    # simple method to coerce the molecule into a neutral charge state.
    # It does this by removing any charge from each atom and then removing the
    # hydrogens (which will then can be added back by the user and will be
    # added back with proper valence).  If the molecule had hydrogens added it
    # will return the molecule with hydrogens added
    # returns self.
    def neutral!
      had_hydrogens = h_added?
      atoms.each {|atom| atom.charge = 0 if (atom.charge != 0) }
      remove_h!
      add_h! if had_hydrogens 
      self
    end

    # adds hydrogens 
    #def add_h_at_ph!(ph=7.4)
    #  # creates a new molecule (currently writes to smiles and uses babel
    #  # commandline to get hydrogens at a given pH; this is because no pH model
    #  # in ruby bindings currently).
    #
    #      ## write the file with the molecule
    #      #self.write_file("tmp.smi")
    #      #system "#{Rubabel::CMD[:babel]} -i smi tmp.smi -p #{ph} -o can tmp.can"
    #      #Molecule.from_file("tmp.can")
    #    end
    #    #alias_method :add_h!, :add_hydrogens!

    # returns self
    def delete_hydrogens!
      @ob.delete_hydrogens
      self
    end
    alias_method :remove_h!, :delete_hydrogens!
    alias_method :remove_hydrogens!, :delete_hydrogens!

    # returns just the smiles string :smi (not the id)
    def smiles
      to_s(:smi)
    end

    # returns just the smiles string (not the id)
    def csmiles
      to_s(:can)
    end

    # checks to see if the molecules are the same OBMol object underneath by
    # modifying one and seeing if the other changes.  This is because
    # openbabel routinely creates new objects that point to the same
    # underlying data store, so even checking for OBMol equivalency is not
    # enough.
    def equal?(other)
      return false unless other.is_a?(self.class)
      are_identical = false
      if self.title == other.title
        begin
          obj_id = self.object_id.to_s
          self.title += obj_id
          are_identical = (self.title == other.title)
        ensure
          self.title.sub(/#{obj_id}$/,'')
        end
        are_identical
      else
        false
      end
    end

    alias_method :eql?, :equal?

    # defined as whether the csmiles strings are identical.  This incorporates
    # more information than the FP2 fingerprint, for instance (try changing
    # the charge and see how it does not influence the fingerprint).
    # Obviously, things like title or data will not be evaluated with ==.  See
    # equal? if you are looking for identity.  More stringent comparisons will
    # have to be done by hand!
    def ==(other)
      other.respond_to?(:csmiles) && (csmiles == other.csmiles)
    end

    # iterates over the molecule's Rubabel::Atom objects
    def each_atom(&block)
      # could use the C++ iterator in the future
      block or return enum_for(__method__)
      iter = @ob.begin_atoms
      atom = @ob.begin_atom(iter)
      while atom
        block.call atom.upcast
        atom = @ob.next_atom(iter)
      end
    end
    alias_method :each, :each_atom

    # iterates over the molecule's Rubabel::Bond objects
    def each_bond(&block)
      # could use the C++ iterator in the future
      block or return enum_for(__method__)
      iter = @ob.begin_bonds
      obbond = @ob.begin_bond(iter)
      while obbond
        block.call obbond.upcast
        obbond = @ob.next_bond(iter)
      end
      self
    end

    # creates a deep copy of the molecule (even the atoms are duplicated)
    def initialize_copy(source)
      super
      @ob = OpenBabel::OBMol.new(source.ob)
      self
    end

    # gets the bond by id
    def bond(id)
      @ob.get_bond_by_id(id).upcast
    end

    # returns the array of bonds.  Consider using #each_bond
    def bonds
      each_bond.map.to_a
    end

    # gets the atom by id
    def atom(id)
      @ob.get_atom_by_id(id).upcast
    end

    # returns the array of atoms. Consider using #each
    def atoms
      each_atom.map.to_a
    end

    def dim
      @ob.get_dimension
    end

    #def each_bond(&block)

    # TODO: implement
    #def calc_desc(descnames=[])
    #end

    # TODO: implement
    # list of supported descriptors. (Not Yet Implemented!)
    #def descs
    #end

    def tanimoto(other, type=DEFAULT_FINGERPRINT)
      other.nil? ? 0 : Rubabel::Molecule.tanimoto(self, other, type)
    end

    # returns a  std::vector<unsigned int> that can be passed directly into
    # the OBFingerprint.tanimoto method
    def ob_fingerprint(type=DEFAULT_FINGERPRINT)
      fprinter = OpenBabel::OBFingerprint.find_fingerprint(type) || raise(ArgumentError, "fingerprint type not found")
      fp = OpenBabel::VectorUnsignedInt.new
      fprinter.get_fingerprint(@ob, fp) || raise("failed to get fingerprint for #{mol}")
      fp
    end

    # obj is an atom or bond
    def delete(obj)
      case obj
      when Rubabel::Bond
        delete_bond(obj)
      when Rubabel::Atom
        delete_atom(obj)
      else 
        raise(ArgumentError, "don't know how to delete objects of type: #{obj.class}")
      end
    end

    # if given a bond, deletes it (doesn't garbage collect).  If given two
    # atoms, deletes the bond between them.
    def delete_bond(*args)
      case args.size
      when 1
        @ob.delete_bond(args[0].ob, false)
      when 2
        @ob.delete_bond(args[0].get_bond(args[1]).ob, false)
      end
    end

    def delete_atom(atom)
      @ob.delete_atom(atom.ob, false)
    end

    # swaps to_move1 for to_move2 on the respective anchors
    # returns self
    def swap!(anchor1, to_move1, anchor2, to_move2)
      OpenBabel::OBBuilder.swap(@ob, *[anchor1, to_move1, anchor2, to_move2].map {|at| at.ob.get_idx } )
      self
    end

    # creates a new (as yet unspecified) bond associated with the molecule and gives it a unique id
    def new_bond
      @ob.new_bond.upcast
    end

    # takes a pair of Rubabel::Atom objects and adds a bond to the molecule
    # returns whether the bond creation was successful.
    def add_bond!(atom1, atom2, order=1)
      @ob.add_bond(atom1.idx, atom2.idx, order)
    end

    # yields self after deleting the specified bonds.  When the block is
    # closed the bonds are restored.  Returns whatever is returned from the
    # block.
    def delete_and_restore_bonds(*bonds, &block)
      bonds.each do |bond|
        unless @ob.delete_bond(bond.ob, false)
          raise "#{bond.inspect} not deleted!" 
        end
      end
      reply = block.call(self)
      bonds.each {|bond| @ob.add_bond(bond.ob) }
      reply
    end

    # splits the molecules at the given bonds and returns the fragments.  Does
    # not alter the caller.  If the molecule is already fragmented, then
    # returns the separate fragments.
    # ##!! Doesn't handle adducts, retained for backwards compatibility
    def basic_split(*bonds)
      if bonds.size > 0
        delete_and_restore_bonds(*bonds) do |mol|
          mol.ob.separate.map(&:upcast)
        end
      else
        self.ob.separate.map(&:upcast)
      end
    end
    # ARRGH!
    def split(*bonds)
      find_internal_adducts(self) if self.adducts.empty?
      to_return = nil
      if bonds.size > 0
        to_return = delete_and_restore_bonds(*bonds) do |mol|
          mol.ob.separate.map(&:upcast)
        end
      else
        to_return = self.ob.separate.map(&:upcast)
      end
      if adducts?
        to_return.map {|a| a.adducts = @adducts}
      end
      to_return
    end

    # Splits, but also provides combined 
    def split_with_adducts(*bonds)
      to_return = split(*bonds)
      if adducts?
        products = to_return.product(@adducts)#.map {|mol, adduct| mol.associate_atom! adduct.atoms.first ; mol}
        products.map! {|mol, adduct| mol2 = mol.dup; mol2.associate_atom! adduct.atoms.first ; mol2}
        to_return = to_return + products
      end
      to_return
    end


    # splits the molecules at the given bonds and returns the fragments.  Does
    # not alter the caller.  If the molecule is already fragmented, then
    # returns the separate fragments.
    def failed_split(*bonds)
      puts "Ground Truth: #{basic_split(*bonds).inspect}"
      self.find_internal_adducts unless @adducts
     # binding.pry if bonds.size == 2
      returns = []
      mols = []
      mols2 = []
      product_mols = []
      if bonds.size > 0
        # set apart the adducts
        p delete_and_restore_bonds(*bonds) do |mol|
          mol.ob.separate.map(&:upcast)
        end
        bonds.each do |bond|
          delete_and_restore_bonds(bond) do |mol| 
            mols = mol.ob.separate.map(&:upcast)
            puts "MOLS: #{mols.inspect}"
            mols = []
            adducts_present = []
            mol.ob.separate.map(&:upcast).map {|a| @adducts.include?(a) ? adducts_present << a : mols << a}
            puts "MOLS: #{mols.inspect}"
            product_mols << mols
          end
        end # bonds.each
        puts "PMOLS: #{product_mols.inspect}"
        mols
      else
        mols = self.ob.separate.map(&:upcast).delete_if {|a| @adducts.include?(a)}
        mols2 = mols.map(&:dup)
        mols_all = mols.map(&:dup)
        @adducts.map do |adduct|
          #mols2.each {|mol| mol.associate_atom! adduct.atoms.first }
          mols_all.product([adduct]).map {|mol, adduct| mol.associate_atom! adduct.atoms.first }
          mols2.product([adduct]).map {|mol, adduct| mol.associate_atom! adduct.atoms.first }
        end
#        p mols3
#        p mols_all
#        puts "+"*50
#        p mols2 != mols
        products = mols2 != mols ? mols.product(mols2) : [mols]
       # p products
        products.delete_if do |set|
          puts "set: #{set}"
          set.last.ob.separate.map(&:upcast).include?(set.first)
          puts "set: #{set}"
        end
       # p products
        #binding.pry
        #puts "adduct_added: #{adduct_added.inspect}"
        # Right now, I'm making the response sets of matched pairs, even if they have adducts... which they should?
        # Is there a better way to feed these back?
        #adduct_added.empty? ? mols : mols.flatten.sort_by(&:mol_wt).reverse.zip(adduct_added.sort_by(&:mol_wt))
        #products.first.is_a?(Array) ? products.flatten : products
        products
      end
    end

    def each_fragment(&block)
      block or return enum_for(__method__)
      @ob.separate.each do |ob_mol|
        block.call( ob_mol.upcast )
      end
    end

    # emits smiles without the trailing tab, newline, or id.  Use write_string
    # to get the default OpenBabel behavior (ie., tabs and newlines).
    def to_s(type=DEFAULT_OUT_TYPE)
      string = write_string(type)
      case type
      when :smi, :smiles, :can
        # remove name with openbabel options in the future
        string.split(/\s+/).first
      else
        string
      end
    end

    # returns a Rubabel::MoleculeData hash
    def data
      Rubabel::MoleculeData.new(@ob)
    end

    # sensitive to add_h!
    def num_atoms(count_implied_hydrogens=false) 
      if !count_implied_hydrogens
        @ob.num_atoms  
      else
        @ob.num_atoms + reduce(0) {|cnt, atom| cnt + atom.ob.implicit_hydrogen_count }
      end
    end
    def num_bonds() @ob.num_bonds  end
    def num_hvy_atoms() @ob.num_hvy_atoms  end
    def num_residues() @ob.num_residues  end
    def num_rotors() @ob.num_rotors  end

    # If filename_or_type is a symbol, then it will return a string of that type.
    # If filename_or_type is a string, will write to the filename
    # given no args, returns a DEFAULT_OUT_TYPE string
    def write(filename_or_type=:can, out_options={})
      if filename_or_type.is_a?(Symbol)
        write_string(filename_or_type, out_options)
      else
        write_file(filename_or_type, out_options)
      end
    end

    # adds hydrogens if necessary.  Performs only steepest descent
    # optimization (no rotors optimized)
    # returns self
    def local_optimize!(forcefield=DEFAULT_FORCEFIELD, steps=500)
      add_h! unless hydrogens_added?
      if dim == 3
        ff = Rubabel.force_field(forcefield.to_s)
        ff.setup(@ob) || raise(OpenBabelUnableToSetupForceFieldError)
        ff.steepest_descent(steps)  # is the default termination count 1.0e-4 (used in obgen?)
        ff.update_coordinates(@ob)
      else
        make_3d!(forcefield, steps) 
      end
      self
    end

    #def global_optimize!(forcefield=DEFAULT_FORCEFIELD, steps=1000)
    #  if dim != 3
    #    # don't bother optimizing yet (steps=nil)
    #    make_3d!(DEFAULT_FORCEFIELD, nil)
    #  end
    #end

    # does a bit of basic local optimization unless steps is set to nil
    # returns self
    def make_3d!(forcefield=DEFAULT_FORCEFIELD, steps=50)
      BUILDER.build(@ob)
      @ob.add_hydrogens(false, true) unless hydrogens_added?
      local_optimize!(forcefield, steps) if steps
      self
    end
    alias_method :make3d!, :make_3d!

    # yields the type of object.  Expects the block to yield the image string.
    def png_transformer(type_s, out_options={}, &block)
      orig_out_options = out_options[:size]
      if type_s == 'png'
        png_output = true
        type_s = 'svg'
        if out_options[:size]
          unless out_options[:size].to_s =~ /x/i
            out_options[:size] = out_options[:size].to_s + 'x' + out_options[:size].to_s
          end
        end
      else
        if out_options[:size].is_a?(String) && (out_options[:size] =~ /x/i)
          warn 'can only use the width dimension for this format'
          out_options[:size] = out_options[:size].split(/x/i).first
        end
      end
      image_blob = block.call(type_s, out_options)
      if png_output
        st = StringIO.new
        image = MiniMagick::Image.read(image_blob, 'svg')
        image.format('png')
        # would like to resize as an svg, then output the png of proper
        # granularity...
        image.resize(out_options[:size]) if out_options[:size]
        image_blob = image.write(st).string
      end
      out_options[:size] = orig_out_options
      image_blob
    end

    # out_options include any of those defined 
    def write_string(type=DEFAULT_OUT_TYPE, out_options={})
      png_transformer(type.to_s, out_options) do |type_s, _out_opts|
        obconv = out_options[:obconv] || OpenBabel::OBConversion.new
        obconv.set_out_format(type_s)
        obconv.add_opts!(:out, _out_opts)
        obconv.write_string(@ob)
      end
    end

    # writes to the file based on the extension given (must be recognized by
    # OpenBabel).  If png is the extension or format, the png is generated
    # from an svg.
    def write_file(filename, out_options={})
      type = Rubabel.filetype(filename)
      File.write(filename, write_string(type, out_options))
    end

    def inspect
      "#<Mol #{to_s}>"
    end

    # returns self
    def convert_dative_bonds!
      @ob.convert_dative_bonds
      self
    end

    # centers the molecule (deals with the atomic coordinate systems for 2D or
    # 3D molecules). returns self.
    def center!
      @ob.center
      self
    end

    # returns self
    def kekulize!
      @ob.kekulize
      self
    end

    # returns self
    def strip_salts!
      @ob.strip_salts!
      self
    end

    # returns self
    def highlight_substructure!(substructure, color='red')
      tmpconv = OpenBabel::OBConversion.new
      tmpconv.add_option("s",OpenBabel::OBConversion::GENOPTIONS, "#{substructure} #{color}")
      self.ob.do_transformations(tmpconv.get_options(OpenBabel::OBConversion::GENOPTIONS), tmpconv)
      self
    end

    ## returns self
    #def highlight_substructure!(substructure)
    #  #obabel benzodiazepine.sdf.gz -O out.svg --filter "title=3016" -s "c1ccc2c(c1)C(=NCCN2)c3ccccc3 red" -xu -d
    #  obconv.set_out_format('svg')

    #  obconv.add_option("s",OpenBabel::OBConversion::GENOPTIONS, "#{substructure} red")
    #  obconv.add_option("d",OpenBabel::OBConversion::GENOPTIONS) 	
    #  self.ob.do_transformations(obconv.get_options(OpenBabel::OBConversion::GENOPTIONS), obconv)

    #  obconv.add_option("u",OpenBabel::OBConversion::OUTOPTIONS)
    #  self
    #end

    def graph_diameter
      distance_matrix = Array.new
      self.atoms.each do |a|
        iter = OpenBabel::OBMolAtomBFSIter.new(self.ob, a.idx)
        while iter.inc.deref do
          distance_matrix << iter.current_depth - 1
        end
      end
      distance_matrix.max
    end
    ADDUCTS_LEGEND.each_pair {|name,str| ADDUCTS[name] = from_string(str)}

  end
end

# NOTES TO SELF:
=begin
~/tools/openbabel-rvmruby1.9.3/bin/babel -L descriptors
abonds    Number of aromatic bonds
atoms    Number of atoms
bonds    Number of bonds
cansmi    Canonical SMILES
cansmiNS    Canonical SMILES without isotopes or stereo
dbonds    Number of double bonds
formula    Chemical formula
HBA1    Number of Hydrogen Bond Acceptors 1 (JoelLib)
HBA2    Number of Hydrogen Bond Acceptors 2 (JoelLib)
HBD    Number of Hydrogen Bond Donors (JoelLib)
InChI    IUPAC InChI identifier
InChIKey    InChIKey
L5    Lipinski Rule of Five
logP    octanol/water partition coefficient
MR    molar refractivity
MW    Molecular Weight filter
nF    Number of Fluorine Atoms
s    SMARTS filter
sbonds    Number of single bonds
smarts    SMARTS filter
tbonds    Number of triple bonds
title    For comparing a molecule's title
TPSA    topological polar surface area

# why can't I find any descriptors????
>> require 'openbabel'
=> true
>> OpenBabel::OBDescriptor.find_type("logP")
=> nil
>> OpenBabel::OBDescriptor.find_type("MR")
=> nil

=end


=begin
  e.g., -xu
 u no element-specific atom coloring
    Use this option to produce a black and white diagram
 U do not use internally-specified color
    e.g. atom color read from cml or generated by internal code
 b black background
    The default is white. The atom colors work with both.
 C do not draw terminal C (and attached H) explicitly
    The default is to draw all hetero atoms and terminal C explicitly,
    together with their attched hydrogens.
 a draw all carbon atoms
    So propane would display as H3C-CH2-CH3
 d do not display molecule name
 e embed molecule as CML
    OpenBabel can read the resulting svg file as a cml file.
 p# scale to bondlength in pixels(single mol only)
 i add index to each atom
    These indices are those in sd or mol files and correspond to the
    order of atoms in a SMILES string.
 j do not embed javascript
    Javascript is not usually embedded if there is one one molecule,
    but it is if the rows and columns have been specified as 1: ``-xr1 -xc1``
 w generate wedge/hash bonds(experimental)
 x omit XML declaration (not displayed in GUI)
    Useful if the output is to be embedded in another xml file.
 A display aliases, if present
    This applies to structures which have an alternative, usually
    shorter, representation already present. This might have been input
    from an A or S superatom entry in an sd or mol file, or can be
    generated using the --genalias option. For example::

      echo "c1cc(C=O)ccc1C(=O)O" | babel -ismi out.svg --genalias -xA

    would add a aliases COOH and CHO to represent the carboxyl and
    aldehyde groups and would display them as such in the svg diagram.
    The aliases which are recognized are in data/superatom.txt, which
    can be edited.

## SPECIFIC TO JUST MULTIPLE MOLS:

 c# number of columns in table
 r# number of rows in table
 N# max number objects to be output
 l draw grid lines

=end

# things JTP is playing with:
=begin
      opts = DEFAULT_DRAW_OPTS.merge( opts )
      self.title = opts[:title] if opts[:title]
      self.obconv.set_out_format(opts[:format].to_s)
      p OpenBabel::OBConversion::OUTOPTIONS
      p OpenBabel::OBConversion::OUTOPTIONS.class
      p self.obconv.get_options(OpenBabel::OBConversion::OUTOPTIONS).methods - Object.new.methods
      abort 'here'
      p self.obconv.get_options(OpenBabel::OBConversion::GENOPTIONS)
      self.obconv.add_option("P", OpenBabel::OBConversion::OUTOPTIONS, opts[:size].to_s) #sets image size to 300X300
      p self.obconv.get_options(OpenBabel::OBConversion::OUTOPTIONS)
      p self.obconv.get_options(OpenBabel::OBConversion::GENOPTIONS)
      p OpenBabel::OBConversion::OUTOPTIONS
      p OpenBabel::OBConversion::OUTOPTIONS.class
      #self.obconv.add_option("gen2D",OpenBabel::OBConversion::GENOPTIONS)
      #self.obconv.add_option("h",OpenBabel::OBConversion::GENOPTIONS)
      #self.ob.do_transformations(self.obconv.get_options(OpenBabel::OBConversion::GENOPTIONS), self.obconv)
      self.obconv.write_file(self.ob, opts[:filename])
=end

#opts = DEFAULT_DRAW_OPTS.merge( opts )
#self.title = opts[:title] if opts[:title]
##self.title = (opts[:title]==nil ? "" : opts[:title])
#self.obconv.set_in_and_out_formats('smi',opts[:format])
#self.obconv.add_option("P",OpenBabel::OBConversion::OUTOPTIONS, opts[:size].to_s) #sets image size to 300X300
#self.obconv.do_transformations(self.obconv
#self.obconv.write_file(self.ob, opts[:filename])

