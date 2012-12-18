%%%Update state
%%%To update a state we need to update the following
%%%X, V1, V2, E1 and dEdE1
function orig_state = update_state(orig_state,update_state,indices, flip_on_rej)
    

orig_state.X(:,indices) = update_state.X(:,indices);
orig_state.V1(:,indices) = update_state.V1(:,indices);
if flip_on_rej == 2
    orig_state.V2(:,indices) = update_state.V2(:,indices);
end
orig_state.E(:,indices) = update_state.E(:,indices);
orig_state.dEdX(:,indices) = update_state.dEdX(:,indices);

end