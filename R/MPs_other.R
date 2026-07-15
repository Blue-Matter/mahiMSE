# Other MPs

# --- Effort at 2022 levels -------------------------------------------
Eff = function(Data, rel_TAE = 1){
    season = get_season(Data)                                          # Get season (this time step, not last with data) [14]
    ad = new('advice')                                                 # Make new Advice object
    ad = do_TAE(Data, TAE, rel_TAE, season, ad)
    ad
}
class(Eff) = "mp"

