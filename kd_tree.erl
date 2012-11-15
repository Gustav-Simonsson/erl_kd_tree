%%%-------------------------------------------------------------------
%%% @author gustav <gustav.simonsson@gmail.com>
%%% @copyright (C) 2012, gustav
%%% @doc
%%% http://en.wikipedia.org/wiki/K-d_tree
%%% @end
%%% Created : 20 Oct 2012 by gustav <gustav.simonsson@gmail.com>
%%%-------------------------------------------------------------------
-module(kd_tree).

-include_lib("proper/include/proper.hrl").

-compile(export_all).

-record(n, {p, l, r}).
-type kd_tree() :: #n{}.
-type point() :: tuple().
-type nearest_neighbour() :: {float(), point()}.

%%%===================================================================
%%% API
%%%===================================================================
%% http://en.wikipedia.org/wiki/K-d_tree
%% K is number of dimensions. Each K-dimensional point below is a tuple with an
%% Id as first element. Dimensional coordinates start at element 2.
%% Example: 2D point: {foo_id, 3, 5}. 3D point: {bar_id, 3, 5, 7}.
kdtree([], _Depth, _K) ->
    nil;
kdtree(Points, Depth, K) ->
    Axis = Depth rem K,
    NewDepth = Depth + 1,
    {Median, {LeftPoints, RightPoints}} = split_on_median_by_axis(Points, Axis),
    #n{p = Median,
        l = kdtree(LeftPoints, NewDepth, K),
        r = kdtree(RightPoints, NewDepth, K)}.

split_on_median_by_axis([Point], _Axis) ->
    {Point, {[], []}};
split_on_median_by_axis(Points, Axis) ->
    Pos = Axis + 2,
    HalfPos = (length(Points) div 2),
    AxisSortedPoints = lists:keysort(Pos, Points),
    Median = lists:nth(HalfPos + 1, AxisSortedPoints),
    {Median, lists:split(HalfPos, lists:delete(Median, AxisSortedPoints))}.

-spec nearest_neighbour(Tree :: kd_tree(), SpatialPoint :: point(),
                        NeighbourCount :: pos_integer(),
                        Dimensions :: pos_integer()) ->
                               NearstNeighbours :: [nearest_neighbour()].

nearest_neighbour(nil, _SearchPoint, _NeighbourCount, _Dimensions) ->
    [];
nearest_neighbour(Tree, SearchPoint, NeighbourCount, Dimensions) ->
    Depth = 0,
    CurrentBests = {NeighbourCount, []},
    {_CBR, CBs} = nn(Tree, SearchPoint, CurrentBests, Depth, Dimensions),
    lists:keysort(1,CBs).

%% Step 2:
nn(#n{l = nil, r = nil, p = P}, SP, {CBR,_} = CBS,  _D, K) when CBR > 0 ->
    add_cb(P, CBS, SP, K);
nn(#n{l = nil, r = nil, p = P}, SP, CBS, _D, K) ->
    update_cb(P, CBS, SP, K);
nn(#n{l = L, r = R, p = P}, SP, CBS1, D, K) ->
    Axis = D rem K,
    %% Step 3:
    {Evaluated, {CBR2,_} = CBS2} =
        case {L,R} of
            %% Step 1:
            {_, nil} -> {all,  nn(L, SP, CBS1, D + 1, K)};
            {nil, _} -> {all, nn(R, SP, CBS1, D + 1, K)};
            _ -> case axis_traverse_left(SP, P, Axis) of
                     true ->  {left, nn(L, SP, CBS1, D + 1, K)};
                     false -> {right, nn(R, SP, CBS1, D + 1, K)}
                 end
        end,
    %% Step 3.1:
    CBS3 = case CBR2 of
               0 -> update_cb(P, CBS2, SP, K);
               _ -> add_cb(P, CBS2, SP, K)
           end,
    %% Step 3.2:
    case shorter_than_some_cb_by_axis(P, CBS3, SP, Axis) of
        true -> case Evaluated of
                    left ->  nn(R, SP, CBS3, D + 1, K);
                    right -> nn(L, SP, CBS3, D + 1, K);
                    all -> CBS3
                end;
        false -> CBS3
    end.

add_cb(P, {CBR,CB}, SP, K) ->
    SQD = sqpd(P, SP, K),
    {CBR - 1, [{SQD,P}|CB]}.

update_cb(P, {CBR,CB}, SP, K) ->
    CBFarthestFirst = lists:sort(fun(A,B) -> A > B end, CB),
    CBReplace =
        fun(_E, {replaced, _} = Acc) -> Acc;
           ({SQD,_P} = E, {_, CBs} = Acc) ->
                case (PSQD = sqpd(SP, P, K)) < SQD of
                    true -> {replaced, [{PSQD, P}|lists:delete(E, CBs)]};
                    false -> Acc
                end
        end,
    {_, NewCBs} =
        lists:foldl(CBReplace, {can_replace, CBFarthestFirst}, CBFarthestFirst),
    {CBR, NewCBs}.

shorter_than_some_cb_by_axis(P, {_CBR,CBs}, SP, Axis) ->
    F = fun({SQD,_P}) -> axis_sqpd(SP, P, Axis) =< SQD end,
    lists:any(F, CBs).

%% One function per K dimensions.
sqpd(Point1, Point2, 2) ->
    K1d = element(2, Point1) - element(2, Point2),
    K2d = element(3, Point1) - element(3, Point2),
    (K1d * K1d) + (K2d * K2d).

axis_sqpd(Point1, Point2, Axis) ->
    D = element(Axis + 2, Point1) - element(Axis + 2, Point2),
    D * D.

axis_traverse_left(Point1, Point2, Axis) ->
    element(Axis + 2, Point1) < element(Axis + 2, Point2).

%%%===================================================================
%%% Internal functions
%%%===================================================================
nearest_neighbour_linear([], _SearchPoint, _Count, _Dimensions) ->
    [];
nearest_neighbour_linear(Coords, SearchPoint, Count, _Dimensions) ->
    SQD = fun(P) -> {d(SearchPoint,P),P} end,
    lists:sublist(lists:keysort(1,lists:map(SQD, Coords)), Count).

d({_Id,X,Y}, {_Id2,X2,Y2}) ->
    Xd = X - X2,
    Yd = Y - Y2,
    Xd*Xd + Yd*Yd.

%%%===================================================================
%%% Tests
%%%===================================================================
proper_test() ->
    proper:quickcheck(?MODULE:prop_nn()).

prop_nn() ->
    Dimensions = 2,
    ?FORALL({X,Y, Coordinates, NeighbourCount},
            {float(), float(), list({dublin, float(), float()}),
             positive_integer()},
            begin
                Tree = kdtree(Coordinates, 0, Dimensions),
                Point = {p,X,Y},
                NNL = nearest_neighbour_linear(Coordinates, Point,
                                               NeighbourCount, Dimensions),
                NNT = nearest_neighbour(Tree, Point,
                                        NeighbourCount, Dimensions),
                lists:sort(NNL) == lists:sort(NNT)
            end).

positive_integer() ->
    ?SUCHTHAT(X, integer(), X > 0).

unit_test({X,Y, Coordinates, NeighbourCount}) ->
    Dimensions = 2,
    Tree = kdtree(Coordinates, 0, Dimensions),
    Point = {p,X,Y},
    NNL = nearest_neighbour_linear(Coordinates, Point,
                                   NeighbourCount, Dimensions),
    NNT = nearest_neighbour(Tree, Point,
                            NeighbourCount, Dimensions),
    io:format("Search point: ~p~n", [Point]),
    io:format("Linear: ~p~n", [NNL]),
    io:format("Tree: ~p~n", [NNT]),
    NNL == NNT.

